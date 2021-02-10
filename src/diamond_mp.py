#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""

import logging
import multiprocessing as mp
import os
import subprocess as sp

# import time
import aux

_logger = logging.getLogger(__name__)


def make_db(input, out):
    # diamond makedb --in <input fasta file> -d <database output>
    try:
        sp.call(['diamond', 'makedb', '--in', input, '-d', out],
                 stderr=sp.DEVNULL, stdin=sp.DEVNULL, stdout=sp.DEVNULL)
    except sp.SubprocessError as e:
        _logger.error('An error occurred while trying to create a DIAMOND database (Input file: {})'.format(input))
        _logger.exception(e)
        return -1
    else:
        return 0


def create_query_(og_list, og_dir, output_dir, query_name):
    query_text = ''
    for og in og_list:
        file = os.path.join(og_dir, og + '.fa')
        query_text += ''.join(aux.change_def_lines(file, og))
        query_text += '\n'

    file_name = query_name + '.fa'
    path = os.path.join(output_dir, file_name)
    with open(path, 'w') as f:
        f.write(query_text)
    return path


def diamond(relax, query, db_in, db_out, output, evalue, create_db, create_query, delete_db, delete_query, query_dir_in,
            query_dir_out, query_file_name=''):

#def diamond(relax, query, db_in, db_out, output, ident, evalue, bitscore, create_db, create_query, delete_db, delete_query, query_dir_in,
#            query_dir_out, query_file_name=''):


    # Check/Create DB
    if create_db:
        db_name = os.path.basename(db_in).split('.')[0]
        _logger.info('Creating DIAMOND database for {}'.format(db_name))
        res = make_db(db_in, db_out)
        if res != 0:
            _logger.warning('The DIAMOND run including {} database did not occur'.format(db_name))
            return -1
    if create_query:
        # in this case the parameter query is a list of Orthogroups and not a path to a fasta file,
        #  so we have to crete the query fasta file and give the path
        n_query = query.split('_')
        if len(n_query) > 1:
            n_query = create_query_(n_query, query_dir_in, query_dir_out, query_file_name)
        else:
            n_query = os.path.join(query_dir_in, query_file_name + '.fa')
            # because in this case we are using the original file
            delete_query = False
    else:
        n_query = query
    # Run DIAMOND
    # diamond blastp -q <query file >—db <database file> —out <output file> —outfmt 6
    db_name = os.path.basename(db_in).split('.')[0]

    _logger.info('Starting DIAMOND search for: {} and {}'.format(db_name, os.path.basename(query).split('.')[0]))

    # Relaxed Search
    if relax:
        sp.call(
            ['diamond', 'blastp', '-q', n_query, '--db', db_out, '--out', output, '--outfmt', '6', 'qseqid', 'sseqid',
             'pident', 'ppos', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
             '--freq-sd', '1000000000', '-e', evalue,
             '--hit-score', '0', '-k', '0'], stderr=sp.DEVNULL, stdin=sp.DEVNULL, stdout=sp.DEVNULL)

#    # Relaxed Search
#    if relax:
#        sp.call(
#            ['diamond', 'blastp', '-q', n_query, '--db', db_out, '--out', output, '--outfmt', '6', 'qseqid', #'sseqid',
#             'pident', 'ppos', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', '--id', #ident, '--query-cover', '30',
#             '--subject-cover', '30', '--min-score', '0', '--freq-sd', '1000000000', '-e', '1000000',
#             '--hit-score', '0', '-k', '0'], stderr=sp.DEVNULL, stdin=sp.DEVNULL, #stdout=sp.DEVNULL)

    # Restrictive search
    else:
        sp.call(
            ['diamond', 'blastp', '-q', n_query, '--db', db_out, '--out', output, '--outfmt', '6', 'qseqid', 'sseqid',
             'pident', 'ppos', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
             '--freq-sd', '1000000000', '-e', evalue,
             '--hit-score', '0', '-k', '1'], stderr=sp.DEVNULL, stdin=sp.DEVNULL, stdout=sp.DEVNULL)

#    # Restrictive search
#    else:
#        sp.call(
#            ['diamond', 'blastp', '-q', n_query, '--db', db_out, '--out', output, '--outfmt', '6', 'qseqid', 'sseqid',
#             'pident', 'ppos', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', '--id', ident, '--query-cover', '40',
#             '--subject-cover', '40', '--min-score', '0', '--freq-sd', '1000000000', '-e', '1000000',
#             '--hit-score', '0', '-k', '1'], stderr=sp.DEVNULL, stdin=sp.DEVNULL, stdout=sp.DEVNULL)

    # Check if the output for this diamond run is available
    try:
        with open(output, 'r'):
            _logger.debug('Output {} opened! '.format(output))
    except IOError as e:
        _logger.warning('Not able to open {}'.format(output))
        _logger.exception(e)
        # Do not return here because we still need to delete the diamond db file

    # Delete database diamond file
    if delete_db:
        try:
            os.remove(db_out + '.dmnd')
        except OSError as e:
            _logger.error('Error trying to delete database {}'.format(db_out + '.dmnd'))
            _logger.exception(e)
    if delete_query:
        try:
            os.remove(query)
        except OSError as e:
            _logger.error('Error trying to delete query {}'.format(query))
            _logger.exception(e)
    return 0


def run(relax, pairs, output, db_storage, n_cpu, evalue, create_db, create_query, delete_db, delete_query,
        query_dir_in='', query_dir_out='', ids_dic={}):

#def run(relax, pairs, output, db_storage, n_cpu, ident, evalue, bitscore, create_db, #create_query, delete_db, delete_query,
#        query_dir_in='', query_dir_out='', ids_dic={}):

    # Get number of cpu to use
    max_cpu = mp.cpu_count()
    # Check if number of given threads is not the higher than all the available ones
    if n_cpu is not None:
        if n_cpu > max_cpu:
            n_cpu = None
    # Define the number of threads to use, if they were not given
    if n_cpu is None:
        if max_cpu > 1:
            n_cpu = max_cpu - 1
        else:
            n_cpu = max_cpu
    _logger.debug('Using {} cpu to DIAMOND runs.'.format(str(n_cpu)))
    for pair in pairs:
        db_name = pair[1].split('/')[-1].split('.')[0]
        db_out = db_storage + db_name
        if isinstance(pair[0], list):
            if len(pair[0]) > 1:
                pair[0] = '_'.join(pair[0])
                query_name = ids_dic[pair[0]]
            else:
                pair[0] = '_'.join(pair[0])
                query_name = pair[0].split('/')[-1].split('.')[0]
        else:
            query_name = pair[0].split('/')[-1].split('.')[0]
        pair.extend([db_out, os.path.join(output, query_name + '|' + db_name), evalue,
                     create_db, create_query, delete_db, delete_query, query_dir_in, query_dir_out, query_name])

#        else:
#            query_name = pair[0].split('/')[-1].split('.')[0]
#        pair.extend([db_out, os.path.join(output, query_name + '|' + db_name), ident, evalue, #bitscore,
#                     create_db, create_query, delete_db, delete_query, query_dir_in, #query_dir_out, query_name])

        # not optimal - should change
        pair.insert(0, relax)

    # Run DIAMOND sequentially in case we have just on pair to run or we don't have enough cores to run in parallel
    if len(pairs) == 1 or n_cpu == 1:
        _logger.debug('Running DIAMOND sequentially.')
        for pair in pairs:
            relax, query, db_in, db_out, output, evalue, create_db, create_query, delete_db, delete_query, query_dir_in, query_dir_out, query_name = pair
            diamond(relax, query, db_in, db_out, output, evalue, create_db, create_query, delete_db, delete_query,
                    query_dir_in, query_dir_out, query_name)

#        for pair in pairs:
#            relax, query, db_in, db_out, output, ident, evalue, bitscore, create_db, create_query, #delete_db, delete_query, query_dir_in, query_dir_out, query_name = pair
#            diamond(relax, query, db_in, db_out, output, ident, evalue, bitscore, create_db, #create_query, delete_db, delete_query,
#                    query_dir_in, query_dir_out, query_name)

    else:
        _logger.debug('Running DIAMOND in parallel.')
        pool = mp.Pool(n_cpu)
        pool.starmap(diamond, pairs)
        pool.close()
        pool.join()

    return 0
