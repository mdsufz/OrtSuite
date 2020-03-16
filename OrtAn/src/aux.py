#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""

import os

def get_reps(file, index):
    res = {}
    with open(file, 'r') as f:
        lines = f.readlines()
    i = 0
    #stop = False
    i_seq = 0
    while i < len(lines): # and not stop:
        if len(index) == 0:
            #stop = True
            break
        if len(lines[i]) > 0:
            if lines[i][0] == '>':
                i_seq += 1
                if i_seq-1 in index:
                    index.remove(i_seq-1)
                    is_seq = True
                    seq = ''
                    def_line = lines[i]
                    i += 1
                    while is_seq and i < len(lines):
                        if len(lines[i]) == 0:
                            i += 1
                        elif lines[i][0] != '>':
                            seq += lines[i]
                            i += 1
                        else:
                            is_seq = False
                    res[def_line] = seq
                else:
                    i += 1
            else:
                i += 1
        else:
            i += 1
    return res


def get_associations_relaxs(dir_in):
    res = {}
    diamond = {}
    files = [f for f in os.listdir(dir_in) if os.path.isfile(os.path.join(dir_in + f)) and not f.startswith('.')]
    for file in files:
        db = file.split('|')[-1]
        with open(os.path.join(dir_in, file)) as f:
            hits = f.readlines()
        for hit in hits:
            # query, target, ident = hit.split()
            query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send = hit.split()
            og = query.split('_')[0]
            if og not in res:
                res[og] = [db]
                diamond[og] = [[db, query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send]]
            elif db not in res[og]:
                res[og].append(db)
                diamond[og].append([db, query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send])
    return res, diamond

def change_def_lines(file, og):
    with open(file, 'r') as f:
        text = f.readlines()
    last = len(text)
    for i in range(last):
        if len(text[i]) > 0 and text[i][0] == '>':
            def_line = text[i].strip()
            new_def_line = def_line[0] + og + '_' + def_line[1:] + '\n'
            text[i] = new_def_line
    return text

def dict_diamond_res(dir, file, multiple_og, previous_dic):
    if not multiple_og:
        og_id, db = file.split('|')
        if og_id not in previous_dic:
            previous_dic[og_id] = []
    with open(os.path.join(dir, file), 'r') as f:
        hits = f.readlines()
    for hit in hits:
        query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send = hit.split()
        if multiple_og:
            og = query.split('_')[0]
            db = file.split('|')[1]
            query = '_'.join(query.split('_')[1:])
            if og not in previous_dic:
                previous_dic[og] = []
            previous_dic[og].append([db, query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send])
        else:
            query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send = hit.split()
            previous_dic[og_id].append([db, query, target, ident, ppos, qlen, slen, qstart, qend, sstart, send])
    return previous_dic


def orthogroups_to_dic(file):
    """
    :param file: Orthogroups.txt file, resulting from OrthoFinder analysis
    :return: dictionary containing the information of the Orthgroups. {OG1: [Seq1, Seq2, ...], OG2 ...}
    """
    res = {}
    single_og = {}
    with open(file, 'r') as f:
        orthogroups = f.readlines()
    for line in orthogroups:
        if len(line) > 0:
            og, prots = line.split(': ', 1)
            res[og.strip()] = [x.strip() for x in prots.split()]
            if len(res[og.strip()]) == 1:
                single_og[og.strip()] = res[og.strip()]
    return res, single_og


def speciesIDs_to_dic(file):
    res = {}
    with open(file, 'r') as f:
        species = f.readlines()
    for line in species:
        if len(line) > 0:
            sID, name = line.split(': ', 1)
            res[sID.strip()] = name.strip()
    return res


def sequencesIDs_to_dic(file):
    res = {}
    with open(file, 'r') as f:
        sequences = f.readlines()
    for line in sequences:
        if len(line) > 0:
            ssIDs, def_line = line.split(': ', 1)
            specieID, sequenceID = ssIDs.split('_')
            specieID, sequenceID = specieID.strip(), sequenceID.strip()
            name = def_line.split(None, 1)[0]
            name = name.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_")
            res[name] = {}
            res[name]['spec'] = specieID
            res[name]['defl'] = def_line
    return res


def overview(total_og, total_kos, og_ko, cog, dog, name):
    with open(name, 'w') as csv_f:

        ko_og = {}
        for og in og_ko:
            for ko in og_ko[og]:
                if ko not in ko_og:
                    ko_og[ko] = []
                ko_og[ko].append(og)

        csv_f.write('Total Orthogroups;' + str(total_og) + '\n')
        csv_f.write('KOs in the database;' + str(total_kos) + '\n')
        # csv_f.write('Relaxed Search \n')
        csv_f.write('(Relaxed Search) Selected orthogroups;' + str(len(og_ko)) + '\n')
        csv_f.write(
            '(Relaxed Search) % Selected orthogroups;' + str(round((len(og_ko) / total_og) * 100, 1)) + '\n')
        csv_f.write('(Relaxed Search) Associated KOs;' + str(len(ko_og)) + '\n')
        csv_f.write('(Relaxed Search) % Associated KOs;' + str(round(((len(ko_og)) / total_kos) * 100, 1)) + '\n')

        kos_per_og = {}
        for og in og_ko:
            n = len(og_ko[og])
            if n not in kos_per_og:
                kos_per_og[n] = 0
            kos_per_og[n] += 1

        og_per_ko = {}
        for ko in ko_og:
            n = len(ko_og[ko])
            if n not in og_per_ko:
                og_per_ko[n] = 0
            og_per_ko[n] += 1

        n_cogs = len(cog)
        n_dogs = len(dog)
        n_dogs_dif = 0
        og_func = {}
        rest_ko = set()
        for og in cog:
            og_func[og] = cog[og]
            rest_ko.update(cog[og])
        for og in dog:
            if len(dog[og]) > 1:
                n_dogs_dif += 1
            og_func[og] = dog[og]
            rest_ko.update(dog[og])

        # csv_f.write('Restrictive Search Results \n')
        csv_f.write('(Restrictive Search) Orthogroups with annotated sequences;' + str(n_cogs + n_dogs) + '\n')
        csv_f.write('(Restrictive Search) % of Orthogroups with annotated sequences;' + str(
            round(((n_cogs + n_dogs) / total_og) * 100, 1)) + '\n')
        csv_f.write('(Restrictive Search) KOs with assigned sequences;' + str(len(rest_ko)) + '\n')
        csv_f.write('(Restrictive Search) % KOs with annotated sequences;' + str(
            round((len(rest_ko) / total_kos) * 100, 1)) + '\n')
        csv_f.write('ConOG;' + str(n_cogs) + '\n')
        csv_f.write('DivOG;' + str(n_dogs) + '\n')
        csv_f.write('DivOG with more than one KO;' + str(n_dogs_dif) + '\n')
        csv_f.write('Relaxed Search to Restrictive Search \n')
        csv_f.write('Lost orthogroups;' + str(len(og_ko) - (n_cogs + n_dogs)) + '\n')
        csv_f.write('% Lost orthogroups;' + str(round(((len(og_ko) - (n_cogs + n_dogs))) / len(og_ko) * 100, 1)) + '\n')
        csv_f.write('Lost KOs;' + str(len(ko_og) - len(rest_ko)) + '\n')
        csv_f.write('% Lost KOs;' + str(round(((len(ko_og) - len(rest_ko)) / len(ko_og)) * 100, 1)) + '\n')
