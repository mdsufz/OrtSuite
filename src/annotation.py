#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import argparse
import json
import logging
import os
import sys

from aux import overview

# from tool import __version__

__author__ = "JoaoSaraiva"
__copyright__ = "JoaoSaraiva"
__license__ = "MIT"

_logger = logging.getLogger(__name__)


def parse_args(args):
    """Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="Runs the annotation step - filter results from restrictive search and return the annotated "
                    "sequences")
    parser.add_argument(
        '-wd',
        '--workingDir',
        dest="WorkingDirectory",
        help="Working Directory",
        required=True
    )
 
#    parser.add_argument(
#        '-ident',
#        '--identity',
#        dest='ident',
#        type=float,
#        help='Identity threshold to filter the diamond results. Default: 40'
#    )
    parser.add_argument(
        '-evalue',
        '--evalue',
        dest='evalue',
        type=float,
        help='Maximum expected value to filter the diamond results. Default: 0.000000001'
    )
    parser.add_argument(
        '-bitscore',
        '--bitscore',
        dest='bitscore',
        type=float,
        help='Required bit-score to filter the diamond results. Default: 50'
    )
#    parser.add_argument(
#        '-qc',
#        '--queryCoverage',
#        dest='query_cov',
#        type=float,
#        help='Query sequence coverage threshold to filter the diamond results. Default: 70'
#    )
#    parser.add_argument(
#        '-sc',
#        '--subjectCoverage',
#        dest='subject_cov',
#        type=float,
#        help='Subject sequence coverage threshold to filter the diamond results. Default: 70'
#    )
#    parser.add_argument(
#        '-ppos',
#        '--percPosMatches',
#        dest='ppos',
#        type=float,
#        help='Percentage of positive matches threshold to filter the diamond results (should be #higher than identity '
#             'threshould). Default: 80'
#)
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
        logfile = os.path.join(args.WorkingDirectory, 'LOGS/annotation.log')
    else:
        logfile = None
    setup_logging(args.loglevel, logfile)

    if not info['restrictive_search']:
        _logger.warning('The Restrictive Search step was not performed correctly.')
        return -1

    # Get the thresholds to use
    _logger.debug('Setting the thresholds to use.')
    """if args.score:
        score_t = args.score
    else:
        score_t = 90.0
    _logger.debug('Score set to {}.'.format(str(score_t)))"""
#    if args.ident:
#        ident_t = args.ident
#    else:
#        ident_t = 40.0
    # if ident_t < float(info['ident_rest']):
    #     ident_t = float(info['ident_rest'])
#    _logger.debug('Identity set to {}.'.format(str(ident_t)))
    if args.evalue:
        evalue_t = args.evalue
    else:
        evalue_t = 0.000000001
    _logger.debug('E-value set to {}.'.format(str(evalue_t)))
    if args.bitscore:
        bitscore_t = args.bitscore
    else:
        bitscore_t = 50.0
    _logger.debug('Bit-score set to {}.'.format(str(bitscore_t)))
#    if args.query_cov:
#        q_cov_t = args.query_cov
#    else:
#        q_cov_t = 70.0
#    _logger.debug('Query sequence coverage set to {}.'.format(str(q_cov_t)))
#    if args.subject_cov:
#        s_cov_t = args.subject_cov
#    else:
#        s_cov_t = 70.0
#    _logger.debug('Subject sequence coverage set to {}.'.format(str(s_cov_t)))
#    if args.ppos:
#        ppos_t = args.ppos
#    else:
#        ppos_t = 80.0
#    _logger.debug('Percentage of positive matches set to {}.'.format(str(ppos_t)))
    _logger.debug('Open results from previous steps.')

    # Restrictive search results
    with open(os.path.join(info['jsons'], 'diamond_res_rest.json'), 'r') as handle:
        rest_res = json.load(handle)
    # Single og results
    with open(os.path.join(info['jsons'], 'single_og_res.json'), 'r') as handle:
        single_res = json.load(handle)
    # Get orthgroups info
    with open(os.path.join(info['jsons'], 'orthogroups.json'), 'r') as handle:
        orthogroups = json.load(handle)
    with open(os.path.join(info['jsons'], 'sequencesIDs.json'), 'r') as handle:
        sequencesIDs = json.load(handle)
    with open(os.path.join(info['jsons'], 'speciesIDs.json'), 'r') as handle:
        speciesIDs = json.load(handle)

    _logger.debug('Filtering and organization of the results.')

    # Choose best hit for repeated querys (note that querys can only be repeated inside the one OG)

    for og in rest_res:
        querys = {}
        # I need the index of the hit to remove in case is necessary
        # hits = rest_res[og].copy()
        index_to_keep = []
        for i in range(len(rest_res[og])):
            query = rest_res[og][i][1]
            # evalue will be the tiebreaker
            evalue = float(rest_res[og][i][3])
            if query not in querys:
                querys[query] = [i, evalue]
                index_to_keep.append(i)
            # in case the current hit has a better value, change the hit
            elif querys[query][1] < evalue:
                # del hits[querys[query][0]]
                index_to_keep.remove(querys[query][0])
                index_to_keep.append(i)
                querys[query] = [i, evalue]

####### CHECK IF WE NEED TO ADD to evalue and bitscore. 
#            # identity will be the tiebreaker
#            ident = float(rest_res[og][i][3])
#            if query not in querys:
#                querys[query] = [i, ident]
#                index_to_keep.append(i)
#            # in case the current hit has a better value, change the hit
#            elif querys[query][1] < ident:
#                # del hits[querys[query][0]]
#                index_to_keep.remove(querys[query][0])
#                index_to_keep.append(i)
#                querys[query] = [i, ident]

        keep_hits = []
        for i in index_to_keep:
            keep_hits.append(rest_res[og][i])
        rest_res[og] = keep_hits

    prot_func = {}
    func_prot = {}
    func_og = {}
    og_prot_func = {}
    og_func_prot = {}
    for og in rest_res:
        og_prot_func[og] = {}
        og_func_prot[og] = {}
        for hit in rest_res[og]:
            db, query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send, evalue, bitscore = hit
            prot_func, func_prot, func_og = filter(db, og, query, target, 
                                                   float(evalue),  float(bitscore), prot_func, func_prot, func_og,
                                                   evalue_t, bitscore_t)

#            db, query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send, evalue, bitscore #= hit
#            prot_func, func_prot, func_og = filter(db, og, query, target, float(ident), float(ppos), #int(qlen),
#                                                   int(slen), int(qstart), int(qend),
#                                                   int(sstart), int(send), float(evalue),  int(bitscore), prot_func, #func_prot, func_og,
#                                                   ident_t, ppos_t, q_cov_t, s_cov_t, evalue_t, bitscore_t)

    # Add results of single og
    for og in single_res:
        og_prot_func[og] = {}
        og_func_prot[og] = {}
        #print(single_res[og])
        for hit in single_res[og]:
            db, q, target, ident, ppos, qlen, slen, qstart, qend, sstart, send, evalue, bitscore = hit
            query = orthogroups[og][0]
            prot_func, func_prot, func_og = filter(db, og, query, target, 
                                                   float(evalue), float(bitscore), prot_func, func_prot, func_og,
                                                   evalue_t, bitscore_t)


#        #print(single_res[og])
#        for hit in single_res[og]:
#            db, q, target, ident, ppos, qlen, slen, qstart, qend, sstart, send, evalue, bitscore = hit
#            query = orthogroups[og][0]
#            prot_func, func_prot, func_og = filter(db, og, query, target, float(ident), float(ppos), #int(qlen),
#                                                   int(slen), int(qstart), int(qend),
#                                                   int(sstart), int(send), float(evalue), int(bitscore), prot_func, #func_prot, func_og,
#                                                   ident_t, ppos_t, q_cov_t, s_cov_t, evalue_t, bitscore_t)

    # add total number of proteins to og and get og_func_prot to use on create_db
    og_to_del = []
    for og in og_prot_func:
        og_prot_func[og]['Total'] = len(orthogroups[og])
        og_prot_func[og]['Unassigned'] = og_prot_func[og]['Total']
        for prot in orthogroups[og]:
            if prot in prot_func:
                for func in prot_func[prot]:
                    if func not in og_prot_func[og]:
                        og_prot_func[og][func] = 0
                        og_func_prot[og][func] = []
                    og_prot_func[og][func] += 1
                    og_prot_func[og]['Unassigned'] -= 1
                    og_func_prot[og][func].append(prot)
        # Remove og with no annotated sequences
        if og_prot_func[og]['Total'] == og_prot_func[og]['Unassigned']:
            og_to_del.append(og)

    for og in og_to_del:
        del og_prot_func[og]


    # Get species with func
    func_spec = {}
    for func in func_prot:
        func_spec[func] = set()
        for prot in func_prot[func]:
            func_spec[func].add(sequencesIDs[prot]['spec'])
        func_spec[func] = list(func_spec[func])

    # List sets to enable json
    for f in func_prot:
        func_prot[f] = list(func_prot[f])

    # Store info to use on the next step
    with open(os.path.join(info['jsons'], 'annotation_og_func_prot.json'), 'w') as handle:
        json.dump(og_func_prot, handle)

    # Write outputs txt
    with open(os.path.join(info['results'], 'Annotation_Protein_Function.txt'), 'w') as f:
        for prot in prot_func:
            for func in prot_func[prot]:
                f.write(prot + '\t' + func + '\n')

    with open(os.path.join(info['results'], 'Annotation_Function_Protein.txt'), 'w') as f:
        for func in func_prot:
            for prot in func_prot[func]:
                f.write(func + '\t' + prot + '\n')

    # Write table of functions assigned to OGs
    functions = list(func_prot.keys())
    with open(os.path.join(info['results'], 'Orthogroups_Annotation.csv'), 'w') as f:
        f.write('OG ID;Total;Unassigned')
        for function in functions:
            f.write(';' + function)

        for og in og_prot_func:
            f.write('\n' + og + ';' + str(og_prot_func[og]['Total']) + ';' + str(og_prot_func[og]['Unassigned']))
            for function in functions:
                if function not in og_prot_func[og]:
                    f.write(';0')
                else:
                    f.write(';' + str(og_prot_func[og][function]))

    # Write table of functions assigned to species
    species = list(speciesIDs.keys())
    with open(os.path.join(info['results'], 'Species_Annotation.csv'), 'w') as f:
        f.write('Specie')
        for specID in species:
            f.write(';' + speciesIDs[specID])

        for func in func_spec:
            f.write('\n' + func)
            for specID in species:
                if specID not in func_spec[func]:
                    f.write(';0')
                else:
                    f.write(';1')

    # Get orthogroups with all the sequences assigned to the same function
    same = {}
    different = {}
    for og in og_prot_func:
        if len(og_prot_func[og]) == 3 and og_prot_func[og]['Unassigned'] == 0:
            same[og] = list(og_prot_func[og].keys())
            same[og].remove('Unassigned')
            same[og].remove('Total')
        else:
            different[og] = list(og_prot_func[og].keys())
            different[og].remove('Unassigned')
            different[og].remove('Total')

    with open(os.path.join(info['results'], 'ConOG.txt'), 'w') as f:
        for og in same:
            f.write(og + '\t' + same[og][0] + '\n')

    with open(os.path.join(info['results'], 'DivOG.txt'), 'w') as f:
        for og in different:
            f.write(og)
            for func in different[og]:
                f.write('\t' + func)
            f.write('\n')

    # Overview file
    # total og
    total_og = info['total_og']
    # total kos
    total_kos = len(info['db_files'])
    # associations
    with open(os.path.join(info['jsons'], 'associations.json'), 'r') as handle:
        associations = json.load(handle)
    overview(total_og, total_kos, associations, same, different, os.path.join(info['results'], 'Overview.csv'))

    # Update information general info dic
    info['annotation'] = True
    # Storing general info dic updated
    with open(os.path.join(info['jsons'], 'general_info.json'), 'w') as handle:
        json.dump(info, handle)

    # complete info table from db with orthogroups
    try:
        with open(os.path.join(info['db_dir'], 'info_db.json'), 'r') as handle:
            db_info = json.load(handle)

        for ko in db_info:
            # Get associated OG
            if ko in func_og:
                for og in func_og[ko]:
                    for i, line in enumerate(db_info[ko]):
                        db_info[ko][i] = [og] + [line]
            else:
                for i, line in enumerate(db_info[ko]):
                    db_info[ko][i] = [''] + [line]

        with open(os.path.join(info['results'], 'orthogroups_more_info_db.csv'), 'w') as f:
            f.write('OG;KO;KO name;Ec number;Reaction')
            for ko in db_info:
                for line in db_info[ko]:
                    f.write('\n' + ';'.join(line))

    except Exception as e:
        _logger.warning('Failed to create orthogroups_more_info_db.csv file')
        _logger.warning(e)

    _logger.info('Done!')

    return 0


def get_coverage(end, start, lenght):
    return ((end-start)/lenght)*100


def filter(func, og, query, target, evalue, bitscore, prot_func, func_prot, func_og,
           evalue_t, bitscore_t):


#def filter(func, og, query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send, prot_func, func_prot, func_og,
#           ident_t, ppos_t, q_cov_t, s_cov_t):
#           ident_t, evalue_t, bitscore_t, ppos_t, q_cov_t, s_cov_t):
#    q_cov = get_coverage(qend, qstart, qlen)
#    s_cov = get_coverage(send, sstart, slen)

#    #mean = ident + ppos + q_cov + s_cov
#    #mean = mean / 4
#    if ident >= ident_t and ppos >= ppos_t and q_cov >= q_cov_t and s_cov >= s_cov_t:
    if evalue <= evalue_t and bitscore >= bitscore_t:
        if query not in prot_func:
            prot_func[query] = set()
        prot_func[query].add(func)
        if func not in func_prot:
            func_prot[func] = set()
        func_prot[func].add(query)
        if func not in func_og:
            func_og[func] = set()
        func_og[func].add(og)
    return prot_func, func_prot, func_og

def run():
    """Entry point for console_scripts
    """
    main(sys.argv[1:])


if __name__ == '__main__':
    run()
