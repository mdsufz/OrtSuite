#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import argparse
import json
import logging
import math
import os
import shutil
import sys
from random import sample

import aux
import diamond_mp

#from tool import __version__

__author__ = "JoaoSaraiva"
__copyright__ = "JoaoSaraiva"
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
        description="Runs the relaxed search step - all the sequences from the database against representative "
                    "sequences from the Orthofinder Orthgroups.")
    parser.add_argument(
        '-wd',
        '--workingDir',
        dest="WorkingDirectory",
        help="Working Directory",
        required=True
    )
    parser.add_argument(
        '-of',
        '--orthofinder',
        dest="orthofinder",
        help="Path to OrthoFinder results directory. Inside that directory should be the Orthgroups, "
             "Orthogroup_Sequences and WorkingDirectory folders",
        required=True
    )
    parser.add_argument(
        '-ident',
        '--identity',
        dest='ident',
        type=int,
        help='Identity threshold to filter the diamond results. DEFAULT: 30'
    )
    parser.add_argument(
        '-t',
        dest='n_cpu',
        type=int,
        help='Number of threads to use in the parallel processing. By default it uses all the cpu available on the '
             'machine'
    )
    parser.add_argument(
        '-del',
        '--delete',
        dest='delete',
        help='To delete the results stored from diamond (use this option if you don\'t want to spend memory space, '
             'between steps)',
        action='store_const',
        const=True
    )
    parser.add_argument(
        '-l',
        '--logfile',
        dest="logfile",
        help="To send log messages to a file in the output directory",
        action='store_const',
        const=True
    )
    parser.add_argument(
        '-v',
        '--verbose',
        dest="loglevel",
        help="set loglevel to DEBUG",
        action='store_const',
        const=logging.DEBUG)
    args = parser.parse_args(args)

    # Validation of all the arguments given
    error_messages = []
    if args.orthofinder:
        if not os.path.isdir(args.orthofinder):
            error_messages.append('ORTHOFINDER - {} does not exist.'.format(args.orthofinder))
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
    except:
        print('Unable to access the information from this project. Are you sure you gave the correct working directory?')
        return -1

    if args.logfile:
        logfile = os.path.join(args.WorkingDirectory, 'LOGS/relaxed_search.log')
    else:
        logfile = None
    setup_logging(args.loglevel, logfile)

    # Check if diamond is executable
    _logger.info('Checking DIAMOND executable')
    if not shutil.which('diamond'):
        _logger.error('Unable to execute diamond command. Please verify if DIAMOND is correctly installed')
        return -1

    # Check if all the OrthoFinder results information is available
    _logger.info('Checking OrthoFinder information')
    # /Orthogroups/Orthgroups.txt
    if not os.path.isfile(os.path.join(args.orthofinder, 'Orthogroups/Orthogroups.txt')):
        _logger.error("Unable to find Orthogroups.txt file. Path {} doesn't exist.".format(os.path.join(args.orthofinder, 'Orthogroups/Orthogroups.txt')))
        return -1
    # /Orthogroup_sequences/
    if not os.path.isdir(os.path.join(args.orthofinder, 'Orthogroup_Sequences/')):
        _logger.error("Unable to find Orthogroup_Sequences directory. Path {} doesn't exist.".format(os.path.join(args.orthofinder, 'Orthogroup_Sequences/')))
        return -1
    # /WorkingDirectory/SpeciesIDs.txt
    if not os.path.isfile(os.path.join(args.orthofinder, 'WorkingDirectory/SpeciesIDs.txt')):
        _logger.error("Unable to find /WorkingDirectory/SpeciesIDs.txt file. Path {} doesn't exist.".format(os.path.join(args.orthofinder, 'WorkingDirectory/SpeciesIDs.txt')))
        return -1
    # /WorkingDirectory/SequencesIDs.txt
    if not os.path.isfile(os.path.join(args.orthofinder,  'WorkingDirectory/SequenceIDs.txt')):
        _logger.error("Unable to find /WorkingDirectory/SequencesIDs.txt file. Path {} doesn't exist.".format(os.path.join(args.orthofinder, 'WorkingDirectory/SequenceIDs.txt')))
        return -1

    # Store OrthoFinder information in json dics
    orthogroups, single_orthogroups = aux.orthogroups_to_dic(
        os.path.join(args.orthofinder, 'Orthogroups/Orthogroups.txt'))
    info['total_og'] = len(orthogroups) + len(single_orthogroups)
    with open(os.path.join(info['jsons'], 'orthogroups.json'), 'w') as handle:
        json.dump(orthogroups, handle)
    with open(os.path.join(info['jsons'], 'single_orthogroups.json'), 'w') as handle:
        json.dump(single_orthogroups, handle)
    speciesIDs = aux.speciesIDs_to_dic(os.path.join(args.orthofinder, 'WorkingDirectory/SpeciesIDs.txt'))
    with open(os.path.join(info['jsons'], 'speciesIDs.json'), 'w') as handle:
        json.dump(speciesIDs, handle)
    sequencesIDs = aux.sequencesIDs_to_dic(os.path.join(args.orthofinder, 'WorkingDirectory/SequenceIDs.txt'))
    with open(os.path.join(info['jsons'], 'sequencesIDs.json'), 'w') as handle:
        json.dump(sequencesIDs, handle)

    # Select representative sequences from each cluster and create the query file for DIAMOND search
    _logger.info('Selecting representative sequences from each Orthogroup')
    querys = {}
    reps_info = {}
    for og in orthogroups:
        n_seqs = len(orthogroups[og])
        # select one representative for each 10 sequences
        n_reps = math.ceil(n_seqs/10)
        rand = sample(range(n_seqs), n_reps)
        # open orthogroup file to get the aa sequences
        file = os.path.join(args.orthofinder, 'Orthogroup_Sequences/', og + '.fa')
        og_res = aux.get_reps(file, rand)
        reps_info[og] = []
        count = 0
        for k in og_res:
            querys['>' + og + '_' + str(count) + '\n'] = og_res[k]
            count += 1
            p_id = k.strip()[1:]
            reps_info[og].append(p_id)
    # Write fasta file with representatives from each cluster
    with open(os.path.join(info['data'], 'query_reps.fa'), 'w') as f:
        for k in querys:
            f.write(k)
            f.write(querys[k])
    # Store information regarding the selected representatives from each orthogroup
    with open(os.path.join(info['jsons'], 'reps_info.json'), 'w') as handle:
        json.dump(reps_info, handle)

    # DIAMOND
    _logger.info('Starting DIAMOND processes')

    if args.n_cpu:
        cpu = args.n_cpu
    else:
        cpu = None

    # Get the identity threshold to use
    if args.ident:
        ident_t = str(args.ident)
    else:
        ident_t = '30'


    # Prepare pairs and parameters for DIAMOND runs
    pairs = []
    query_file = os.path.join(info['data'], 'query_reps.fa')
    for file in info['db_files']:
        pairs.append([query_file, os.path.join(info['db_dir'], file)])
    output = info['diamond_res_relax']
    db_storage = info['diamond_dbs']

    if args.delete:
        del_db = True
    else:
        del_db = False

    # DIAMOND runs
    diamond_mp.run(True, pairs, output, db_storage, cpu, ident_t,
                   create_db=True, create_query=False, delete_db=del_db, delete_query=False)

    # Get associations between OG and functions
    _logger.info('Storing the associations between the Orthogroups and DB functions')
    # pprint.pprint(aux.get_associations_relaxs(output))
    associations, diamond_res = aux.get_associations_relaxs(output)
    # store json with associations
    with open(os.path.join(info['jsons'], 'associations.json'), 'w') as handle:
        json.dump(associations, handle)
    # create txt with the results
    with open(os.path.join(info['results'], 'Associations_Relaxed_S.txt'), 'w') as f:
        f.write('Orthogroup\tFunction\n')
        for assoc in associations:
            for func in associations[assoc]:
                f.write(assoc + '\t' + func + '\n')

    # Store DIAMOND results for single Orthogroups
    single_og_res = {}
    for og in single_orthogroups:
        if og in diamond_res:
            single_og_res[og] = diamond_res[og]
    with open(os.path.join(info['jsons'], 'single_og_res.json'), 'w') as handle:
        json.dump(single_og_res, handle)

    # Delete DIAMOND result files

    to_delete = [x for x in os.listdir(info['diamond_res_relax']) if os.path.isfile(
        os.path.join(info['diamond_res_relax'], x))]
    for file in to_delete:
        try:
            os.remove(os.path.join(info['diamond_res_relax'], file))
        except OSError as e:
            _logger.error('Unable to delete {}'.format(os.path.join(info['diamond_res_relax'], file)))
            _logger.exception(e)

    # Update info about the project
    info['relaxed_search'] = True
    info['og_sequences_dir'] = os.path.join(args.orthofinder, 'Orthogroup_Sequences/')
    info['need_to_create_db'] = del_db
    info['ident_relax'] = ident_t

    with open(os.path.join(info['jsons'], 'general_info.json'), 'w') as handle:
        json.dump(info, handle)

    _logger.info('Done!')
    return 0

def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])

if __name__ == '__main__':
    run()
