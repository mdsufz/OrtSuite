#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import argparse
import json
import logging
import os
# import subprocess
import shutil
import sys

import aux
import diamond_mp

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
        description="Runs the restrictive search step - all the sequences from the selected Orthgroups versus the "
                    "functions associated with.")
    parser.add_argument(
        '-wd',
        '--workingDir',
        dest="WorkingDirectory",
        help="Working Directory",
        required=True
    )
    parser.add_argument(
        '-t',
        dest='n_cpu',
        type=int,
        help='Number of threads to use in the parallel processing. By default it uses all the cpu available on the '
             'machine.'
    )
    parser.add_argument(
        '-ident',
        '--identity',
        dest='ident',
        type=int,
        help='Identity threshold to filter the diamond results.'
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
        const=logging.DEBUG
    )
    args = parser.parse_args(args)
    if not os.path.isdir(args.WorkingDirectory):
        print('Working Directory - {} does not exist.'.format(args.WorkingDirectory))
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
        print(
            'Unable to access the information from this project. Are you sure you gave the correct working directory?')
        return -1

    if args.logfile:
        logfile = os.path.join(args.WorkingDirectory, 'LOGS/restrictive_search.log')
    else:
        logfile = None
    setup_logging(args.loglevel, logfile)

    if not info['relaxed_search']:
        _logger.warning('The Relaxed Search step was not performed correctly.')
        return -1

    # Check if diamond is executable
    _logger.info('Checking DIAMOND executable')
    if not shutil.which('diamond'):
        _logger.error('Unable to execute diamond command. Please verify if DIAMOND is correctly installed')
        return -1

    db_og = {}
    # with open(os.path.join(info['jsons'], 'orthogroups.json'), 'r') as handle:
    #     orthogroups = json.load(handle)
    with open(os.path.join(info['jsons'], 'single_orthogroups.json'), 'r') as handle:
        single_orthogroups = json.load(handle)
    with open(os.path.join(info['jsons'], 'associations.json'), 'r') as handle:
        associations = json.load(handle)
    # for each DB function get the associated orthogroups
    for og in associations:
        if og not in single_orthogroups:
            for db in associations[og]:
                if db not in db_og:
                    db_og[db] = [og]
                else:
                    db_og[db].append(og)

    if args.n_cpu:
        cpu = args.n_cpu
    else:
        cpu = None

    # Get the identity threshold to use
    if args.ident:
        ident_t = str(args.ident)
    else:
        ident_t = '90'

    if float(ident_t) < float(info['ident_relax']):
        ident_t = info['ident_relax']


    pairs = []
    query_files_ids = {}
    dif_results_files = []
    id_count = 1
    for db in db_og:
        db_file = os.path.join(info['db_dir'], db + '.fa')
        pairs.append([db_og[db], db_file])
        # Create ids to name query files
        if len(db_og[db]) > 1:
            ogs_name = '_'.join(db_og[db])
            if ogs_name not in query_files_ids:
                query_files_ids['_'.join(db_og[db])] = str(id_count)
            dif_results_files.append(query_files_ids[ogs_name] + '|' + db)
            id_count += 1
    diamond_mp.run(False, pairs, info['diamond_res_rest'], info['diamond_dbs'], cpu, ident_t,
                   create_db=info['need_to_create_db'], create_query=True, delete_db=False, delete_query=False,
                   query_dir_in=info['og_sequences_dir'], query_dir_out=info['data'], ids_dic=query_files_ids)
    # diamond_mp.run(pairs, info['diamond_res_rest'], info['diamond_dbs'], cpu, ident_t,
    #               create_db=True, create_query=True, delete_db=True, delete_query=False,
    #               query_dir_in=info['og_sequences_dir'], query_dir_out=info['data'], ids_dic=query_files_ids)

    # Analise results
    res = {}
    for file in os.listdir(info['diamond_res_rest']):
        if file in dif_results_files:
            res = aux.dict_diamond_res(info['diamond_res_rest'], file, True, res)
        else:
            res = aux.dict_diamond_res(info['diamond_res_rest'], file, False, res)
    # Storing results in a json file
    with open(os.path.join(info['jsons'], 'diamond_res_rest.json'), 'w') as handle:
        json.dump(res, handle)

    # Delete DIAMOND result files
    to_delete = [x for x in os.listdir(info['diamond_res_rest']) if os.path.isfile(
        os.path.join(info['diamond_res_rest'], x))]
    for file in to_delete:
        try:
            os.remove(os.path.join(info['diamond_res_rest'], file))
        except OSError as e:
            _logger.error('Unable to delete {}'.format(os.path.join(info['diamond_res_rest'], file)))
            _logger.exception(e)

    # Updating general info dic
    info['restrictive_search'] = True
    info['ident_rest'] = ident_t
    # Storing general info dic updated
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