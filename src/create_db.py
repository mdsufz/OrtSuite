#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import argparse
import json
import logging
import os
import shutil
import sys

# from tool import __version__

__author__ = "MartaLopesGomes"
__copyright__ = "MartaLopesGomes"
__license__ = "mit"

_logger = logging.getLogger(__name__)


def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Runs the create database step - the annotated sequences are added to the previous database.")
    group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument(
        '-wd',
        '--workingDir',
        dest="WorkingDirectory",
        help="Working Directory",
        required=True
    )
    group.add_argument(
        '-o',
        '--out',
        dest="output",
        help="Output directory to create new database (initial database + new annotated sequences)."
    )
    group.add_argument(
        '-up',
        '--update',
        dest="update",
        help="Use this option to update the initial database with the new annotated sequences. "
             "ATTENTION: If you use this option the initial database will be changed permanently",
        action="store_const",
        const=True
    )
    parser.add_argument(
        '-l',
        '--logfile',
        dest="logfile",
        help="To send log messages to a file int the output directory",
        action='store_const',
        const=True
    )
    parser.add_argument(
        '-v',
        '--verbose',
        dest="loglevel",
        help="set loglevel to DEBUG",
        action='store_const',
        const=logging.DEBUG
    )
    args = parser.parse_args(args)
    error_messages = []
    if not os.path.isdir(args.WorkingDirectory):
        error_messages.append('Working Directory - {} does not exist.'.format(args.WorkingDirectory))
    if args.output:
        if not os.path.isdir(args.output):
            error_messages.append('OUTPUT - {} does not exist.'.format(args.output))
    if len(error_messages) > 0:
        print('INPUT ARGUMENTS ERROR:')
        for error in error_messages:
            print(error)
            print()
        parser.print_help()
        return None

    return args


def setup_logging(loglevel, logfile):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    # By default I want the logging to show all the information from level INFO
    if not loglevel:
        loglevel = logging.INFO
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(filename=logfile, level=loglevel, format=logformat, datefmt="%Y-%m-%d %H:%M:%S")


def main(args):
    """Main entry point allowing external calls

    Args:
        args ([str]): command line parameter list
    """
    args = parse_args(args)
    if args is None:
        return -1

    # Check project info
    try:
        with open(os.path.join(args.WorkingDirectory, 'jsons/general_info.json'), 'r') as handle:
            info = json.load(handle)
    except IOError as e:
        print(
            'Unable to access the information from this project. Are you sure you gave the correct working directory?')
        print(e)
        return -1

    if args.logfile:
        logfile = os.path.join(args.WorkingDirectory, 'LOGS/create_db.log')
    else:
        logfile = None
    setup_logging(args.loglevel, logfile)

    if not info['annotation']:
        _logger.warning('The Annotation step was not performed correctly.')
        return -1

    with open(os.path.join(info['jsons'], 'annotation_og_func_prot.json'), 'r') as handle:
        og_func_prot = json.load(handle)

    with open(os.path.join(info['jsons'], 'sequencesIDs.json'), 'r') as handle:
        sequencesIDs = json.load(handle)

    if args.update:
        db_dir = info['db_dir']
    else:
        db_dir = args.output
        for file in info['db_files']:
            try:
                shutil.copy(os.path.join(info['db_dir'], file), os.path.join(db_dir, file))
            except IOError as e:
                _logger.warning('Unable to copy {} from {} to {}.'.format(file, info['db_dir'], db_dir))
                _logger.exception(e)

    _logger.info('Collecting sequences from orthogroups fasta files.')
    # Get sequences to add (from each og)
    sequences = {}
    func_prots = {}
    for og in og_func_prot:
        seqs_to_get = set()
        for func in og_func_prot[og]:
            if func not in func_prots:
                func_prots[func] = set()
            seqs_to_get.update(og_func_prot[og][func])
            func_prots[func].update(og_func_prot[og][func])
        file = os.path.join(info['og_sequences_dir'], og + '.fa')
        _logger.debug('Working on {}'.format(file))
        sequences = get_sequences(sequences, file, seqs_to_get)

    _logger.info('Adding sequences to the database.')
    # Add sequences
    for func in func_prots:
        seqs = {}
        file = os.path.join(db_dir, func + '.fa')
        _logger.debug('Working on {}'.format(file))
        for seqID in func_prots[func]:
            seqs[sequencesIDs[seqID]['defl']] = sequences[seqID]
        add_sequences(file, seqs)

    _logger.info('Writing output.')
    # Create tracking file of the added sequences to each function file
    with open(os.path.join(info['results'], 'db_tracking.txt'), 'w') as f:
        f.write('Function\tTotal seqs added\tsequences')
        for func in func_prots:
            f.write('\n'+func+'\t'+str(len(func_prots[func])))
            for prot in func_prots[func]:
                f.write('\t'+prot)
    _logger.info('Done!')
    return 0


def get_sequences(dic, file, ids):
    with open(file, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines) and len(ids) > 0:
        if len(lines[i]) > 0:
            if lines[i][0] == '>':
                name = lines[i][1:].strip()
                if name in ids:
                    ids.remove(name)
                    is_seq = True
                    seq = ''
                    i += 1
                    while is_seq and i < len(lines):
                        if len(lines[i]) == 0:
                            i += 1
                        elif lines[i][0] != '>':
                            seq += lines[i]
                            i += 1
                        else:
                            is_seq = False
                    dic[name] = seq
                else:
                    i += 1
            else:
                i += 1
        else:
            i += 1
    return dic


def add_sequences(file, sequences):
    text = ''
    for seq in sequences:
        text += '\n' + '>' + seq + '\n' + sequences[seq]
    with open(file, 'a') as f:
        f.write(text)
    return 0


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == '__main__':
    run()