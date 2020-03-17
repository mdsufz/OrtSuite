#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import argparse
import json
import logging
import os
import sys

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
        description="Creates a new project and returns the working directory to use in the other steps of the pipline. "
                    "(Necessary to access the temporary data that the tool stores)")
    # parser.add_argument(
    #    '-gn',
    #    '--genomes',
    #    dest="genomes",
    #    help="Path to folder containing the genomes in separated fasta files",
    #    required=True
    # )
    parser.add_argument(
        '-out',
        '--output',
        dest="output",
        help="Path to folder to store the output results",
        required=True
    )
    parser.add_argument(
        '-db',
        '--data',
        dest="data",
        help="Path to folder containing the fasta files with sequences from the functions of interest",
        required=True
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
    # Validation of all the arguments given
    error_messages = []
    if args.output:
        if not os.path.isdir(args.output):
            message = 'OUTPUT - {} does not exist.'.format(args.output)
            error_messages.append(message)
    if args.data:
        if not os.path.isdir(args.data):
            message = 'DATA - {} does not exist.'.format(args.data)
            error_messages.append(message)
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


# Check if files have fasta extensions and rename them if they have characters :, | or white space
def check_files(folder):
    def is_fasta(f):
        fasta_ext = [".fasta", ".faa", ".fa", ".fas", ".fsa"]
        for ext in fasta_ext:
            if ext in f:
                return True
        return False

    files = [f for f in os.listdir(folder) if is_fasta(f)]
    pro_char = [':', '|', ' ']
    for f in files.copy():
        for char in pro_char:
            if char in f:
                new_name = f.replace(':', '_').replace('|', '_').replace(' ', '_').split('.')[0]
                new_name += '.fa'
                try:
                    os.rename(folder + f, folder + new_name)
                    files.remove(f)
                    files.append(new_name)
                except OSError as e:
                    _logger.warning('Unable to rename a file {} from the database. This could lead to problems in the '
                                    'future steps.'.format(f))
                    _logger.error(e)
                continue
    return files


def main(args):
    """Main entry point allowing external calls

    Args:
      args ([str]): command line parameter list
    """
    args = parse_args(args)
    if args is None:
        return -1
    WorkingDirectory = args.output

    # Working directory structure
    # WorkingDirectory (or the name given from the user)
    #     |_Results
    #     |_jsons
    #     |_LOGS
    #     |_data
    #       |_diamond_dbs
    #       |_diamond_results
    #           |_relax
    #           |_rest

    # Create directories
    try:
        os.mkdir(os.path.join(WorkingDirectory, 'Results/'))
        os.mkdir(os.path.join(WorkingDirectory, 'jsons/'))
        os.mkdir(os.path.join(WorkingDirectory, 'LOGS/'))
        os.mkdir(os.path.join(WorkingDirectory, 'data/'))
        os.mkdir(os.path.join(WorkingDirectory, 'data/diamond_dbs/'))
        os.mkdir(os.path.join(WorkingDirectory, 'data/diamond_results/'))
        os.mkdir(os.path.join(WorkingDirectory, 'data/diamond_results/relax/'))
        os.mkdir(os.path.join(WorkingDirectory, 'data/diamond_results/rest/'))
    except OSError as e:
        print('Error trying to create new directory on : {}'.format(WorkingDirectory))
        print(e)
        return -1

    if args.logfile:
        logfile = os.path.join(WorkingDirectory, 'LOGS/create_project.log')
    else:
        logfile = None
    setup_logging(args.loglevel, logfile)

    # Check input functions
    function_files = check_files(args.data)
    general_info = {'relaxed_search': False,
                    'restrictive_search': False,
                    'annotation': False,
                    'db_dir': args.data,
                    'db_files': function_files,
                    'results': os.path.join(WorkingDirectory, 'Results/'),
                    'jsons': os.path.join(WorkingDirectory, 'jsons/'),
                    'logs': os.path.join(WorkingDirectory, 'LOGS/'),
                    'data': os.path.join(WorkingDirectory, 'data/'),
                    'diamond_res_relax': os.path.join(WorkingDirectory, 'data/diamond_results/relax/'),
                    'diamond_res_rest': os.path.join(WorkingDirectory, 'data/diamond_results/rest/'),
                    'diamond_dbs': os.path.join(WorkingDirectory, 'data/diamond_dbs/')}

    with open(os.path.join(WorkingDirectory, 'jsons/', 'general_info.json'), 'w') as handle:
        json.dump(general_info, handle)
    _logger.info('The results of this project will be stored at: {}'.format(WorkingDirectory))

    return 0


def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])

if __name__ == '__main__':
    run()