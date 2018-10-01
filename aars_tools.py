"""Python tools for aaRS data

a.k.a. "Paul's tools"
"""
import sys
import os.path
from itertools import product
import glob
import json
import xml.etree.ElementTree as ET
import numpy as np
from Bio.Seq import translate

try:
    from levenshtein import levenshtein_distance_c as levenshtein_distance
except ImportError:
    def levenshtein_distance(source, target):
        """Custom Levenshtein distance calculation"""
        distance_matrix = np.zeros((len(source) + 1, len(target) + 1), int)
        distance_matrix[0, :] = np.arange(len(target) + 1)
        distance_matrix[:, 0] = np.arange(len(source) + 1)
        for i, j in product(range(1, len(source)+1), range(1, len(target)+1)):
            substitution_cost = 0 if source[i-1] == target[j-1] else 1
            distance_matrix[i, j] = min(
                distance_matrix[i-1, j] + 1,
                distance_matrix[i, j-1] + 1,
                distance_matrix[i-1, j-1] + substitution_cost
            )
        return int(distance_matrix[len(source), len(target)])


AARS_XML = 'AARS.xml'
BASE_DIR = 'BEAST 2/XMLs/Better Priors (final, actually used XMLs)/'
A_B_E_ZIP = '../../Documents'
CLASS_I = BASE_DIR + 'ClassI_betterPriors.xml'
CLASS_II = BASE_DIR + 'ClassII_betterPriors.xml'


def final_sequences():
    """dictionary of final sequences (no gaps)"""
    tree = ET.parse(AARS_XML)
    root = tree.getroot()
    data = root[0]
    return {s.attrib['taxon']: s.attrib['value']
            for s in data if s.tag == 'sequence'}


def all_sequences():
    """Union of class I and class II sequences"""
    return dict(class_one_sequences(), **class_two_sequences())


def class_one_sequences():
    """dictionary of class I sequences (gaps)"""
    tree = ET.parse(CLASS_I)
    root = tree.getroot()
    data = root[0]
    return {s.attrib['taxon']: s.attrib['value']
            for s in data if s.tag == 'sequence'}


def class_two_sequences():
    """dictionary of class II sequences (gaps)"""
    tree = ET.parse(CLASS_II)
    root = tree.getroot()
    data = root[0]
    return {s.attrib['taxon']: s.attrib['value']
            for s in data if s.tag == 'sequence'}


def find_match(seq_name, no_gaps, gaps):
    """Find best gapped sequences matching the named sequence"""
    try:
        source = no_gaps[seq_name]
    except ValueError:
        print('sequences name', seq_name, 'not found in final sequences file')
        return None
    best_key, best_value = sorted(
        gaps.items(),
        key=(lambda target: levenshtein_distance(source, target[1])))[0]
    return best_key, levenshtein_distance(source, best_value)


def parse_no_gaps_keys(no_gaps_key):
    """Parse key into domain, species, amino acid"""
    domain_dict = {'arch': 'Archaea',
                   'bact': 'Bacteria',
                   'euk': 'Eukaryotes'}
    splits = no_gaps_key.split('_')
    aa = splits[0]
    try:
        domain = domain_dict[splits[1]]
    except KeyError:
        domain = 'Unknown'
    for i, text in enumerate(splits):
        if len(text) == 1 and text[0] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            species = '_'.join(splits[i:])
            break
    else:
        species = splits[-1]
    return domain, species, aa


def make_no_gaps_to_gaps_file():
    """Write the no gaps to gaps file"""
    no_gaps = final_sequences()
    gaps = all_sequences()
    sys.stdout.write('0%')
    sys.stdout.flush()
    data = {}
    for i, no_gaps_key in enumerate(no_gaps.keys()):
        gaps_key, distance = find_match(
            no_gaps_key, no_gaps=no_gaps, gaps=gaps)
        domain, species, aa = parse_no_gaps_keys(no_gaps_key)
        path, path_cost = predict_path(domain, species, aa)
        aa_paths = predict_amino_path(path, aa)
        nuc_paths = predict_nucleotide_path(path, aa)
        aa_data = [read_path_data(p) for p in aa_paths]
        aa_align = [align(gaps[gaps_key], d) for d in aa_data]
        aa_visual = [''.join([c1.lower() if c2 == '-' else c1
                              for c1, c2 in zip(dat, ali)])
                     for dat, ali in zip(aa_data, aa_align)]
        nuc_data = [read_path_data(p) for p in nuc_paths]
        nuc_trans = [translate(p) for p in nuc_data]
        nuc_trans_align = [align(gaps[gaps_key], d) for d in nuc_trans]
        nuc_align = [''.join(['---' if a[i] == '-' else n[i*3:(i+1)*3]
                              for i in range(len(a))])
                     for a, n in zip(nuc_trans_align, nuc_data)]
        nuc_visual = [''.join([n[i*3:(i+1)*3].lower() if a[i] == '-' else n[i*3:(i+1)*3]
                               for i in range(len(a))])
                      for a, n in zip(nuc_trans_align, nuc_data)]
        data[no_gaps_key] = {
            'no-gaps-key': no_gaps_key,
            'domain': domain,
            'species': species,
            'amino acid': aa,
            'gaps-key': gaps_key,
            'no-gaps-value': no_gaps[no_gaps_key],
            'gaps-value': gaps[gaps_key],
            'gaps-sequence-uncertainty': distance / len(no_gaps[no_gaps_key]),
            'base-path': path,
            'path-uncertainty': path_cost / len(species),
            'amino-acid-paths': aa_paths,
            'nucleotide-paths': nuc_paths,
            'amino-acid-data': aa_data,
            'amino-acid-aligned-gaps': aa_align,
            'amino-acid-visual-alignment': aa_visual,
            'nucleotide-data': nuc_data,
            'nucleotide-translation': nuc_trans,
            'nucleotide-trans-align': nuc_trans_align,
            'nucleotide-aligned-gaps': nuc_align,
            'nucleotide-visual-alignment': nuc_visual
        }
        sys.stdout.write('\r{}%'.format((i * 100) // len(no_gaps)))
        sys.stdout.flush()
    with open('gap_data.json', 'w') as f:
        json.dump(data, f, indent=2)


def predict_amino_path(path, aa):
    """Try to find the amino acid sequence file"""
    if path[-3:] == "ile" and aa == 'ile':  # M_modile ends in 'ile' but may not be aa 'ile'
        return predict_amino_path(path, '_ile')
    if 'Musca domestica' in path and aa == 'phe':  # This one has missing resources
        try:  # try using the ALPHA and BETA instead
            paths = (predict_amino_path(path, 'pheALPHA') +
                     predict_amino_path(path, 'pheBETA'))
        except RuntimeError:
            try:
                # phe sometimes has pheALPHA but not pheBETA
                paths = predict_amino_path(path, 'pheALPHA')
            except RuntimeError:
                # phe sometimes has pheBETA but not pheALPHA
                paths = predict_amino_path(path, 'pheBETA')
        return paths
    paths = glob.glob(path + '/amino*/*{}_*'.format(aa))
    if not paths:
        paths = glob.glob(path + '/**/*{}_aa*'.format(aa), recursive=True)
    if not paths:
        if aa == 'glu':  # glu is sometimes listed as gluPro
            return predict_amino_path(path, 'gluPro')
        if aa == 'leu':  # leu is sometimes listed as leuALPHA and leuBETA
            return (predict_amino_path(path, 'leuALPHA') +
                    predict_amino_path(path, 'leuBETA'))
        if aa == 'met':  # met is sometimes listed as leuALPHA and leuBETA
            return (predict_amino_path(path, 'metALPHA') +
                    predict_amino_path(path, 'metBETA'))
        if aa == 'phe':  # phe is sometimes listed as pheALPHA and pheBETA
            try:
                paths = (predict_amino_path(path, 'pheALPHA') +
                         predict_amino_path(path, 'pheBETA'))
            except RuntimeError:
                try:
                    # phe sometimes has pheALPHA but not pheBETA
                    paths = predict_amino_path(path, 'pheALPHA')
                except RuntimeError:
                    # phe sometimes has pheBETA but not pheALPHA
                    paths = predict_amino_path(path, 'pheBETA')
            return paths
        if aa == 'tyr':  # tyr is sometimes listed as Tyr
            return predict_amino_path(path, 'Tyr')
        if aa == 'val':  # val is sometimes listed as Val or val1
            try:
                paths = predict_amino_path(path, 'Val')
            except RuntimeError:
                paths = predict_amino_path(path, 'val1')
            return paths

    if not paths:
        # exceptions
        if "Pyrococcus horikoshii" in path and 'asn' in aa:
            print('\nPyrococcus horikoshii has no asn amino acid sequence')
        elif "Homo sapiens" in path and 'phe' in aa:
            print('\nHomo sapiens has no phe amino acid sequence')
        else:
            raise RuntimeError('Could not find amino acid path for {} ({})'.format(
                os.path.basename(path), aa))
    return paths


def predict_nucleotide_path(path, aa):
    """Try to find the nucleotide sequence file"""
    if path[-3:] == "ile" and aa == 'ile':  # M_modile ends in 'ile' but may not be aa 'ile'
        return predict_nucleotide_path(path, '_ile')
    if 'Methanopyrus kandleri' in path and aa == 'arg':
        return []  # bad data in Mkandleri_arg_nuc
    if 'Crassostrea gigas' in path and aa == 'ala':
        return []  # aa data in Cgigas_ala_nuc (should be nuc)
    if 'Halobacterium sp.' in path and aa == 'gly':
        return []  # bad data format
    if 'Phycisphaera mikurensis' in path and aa == 'lys':
        return []  # aa data in Pmikurensis_lys_nuc (should be nuc)
    paths = glob.glob(path + '/nuc*/*{}_*'.format(aa))
    if not paths:
        paths = glob.glob(path + '/**/*{}_nuc*'.format(aa), recursive=True)
    if not paths:
        if aa == 'glu':  # glu is sometimes listed as gluPro
            return predict_nucleotide_path(path, 'gluPro')
        if aa == 'leu':  # leu is sometimes listed as leuALPHA and leuBETA
            return (predict_nucleotide_path(path, 'leuALPHA') +
                    predict_nucleotide_path(path, 'leuBETA'))
        if aa == 'met':  # met is sometimes listed as leuALPHA and leuBETA
            return (predict_nucleotide_path(path, 'metALPHA') +
                    predict_nucleotide_path(path, 'metBETA'))
        if aa == 'phe':  # phe is sometimes listed as pheALPHA and pheBETA
            try:
                paths = (predict_nucleotide_path(path, 'pheALPHA') +
                         predict_nucleotide_path(path, 'pheBETA'))
            except RuntimeError:
                try:
                    # phe sometimes has pheALPHA but not pheBETA
                    paths = predict_nucleotide_path(path, 'pheALPHA')
                except RuntimeError:
                    # phe sometimes has pheBETA but not pheALPHA
                    paths = predict_nucleotide_path(path, 'pheBETA')
            return paths
    if not paths:
        # exceptions
        if "Bos taurus" in path and 'leu' in aa:
            print('\nBos taurus has no leu nucleotide sequence')
        elif "Pyrococcus horikoshii" in path and 'asn' in aa:
            print('\nPyrococcus horikoshii has no asn nucleotide sequence')
        elif "Giardia lamblia" in path and 'asn' in aa:
            print('\nGiardia lamblia has no asn nucleotide sequence')
        elif "Homo sapiens" in path and 'phe' in aa:
            print('\nHomo sapiens has no phe nucleotide sequence')
        else:
            raise RuntimeError('Could not find nucleotide path for {} ({})'.format(
                os.path.basename(path), aa))
    return paths


def predict_path(domain, species, aa):
    """Try to find the path to the original sequences"""
    base_dir = A_B_E_ZIP + "/Archaeans_Bacteria_Eukaryotes"
    possibles = glob.glob(base_dir + '/{}/{}*'.format(domain, species[0]))
    if not possibles:
        possibles = glob.glob(base_dir + '/{}/*'.format(domain))
    if not possibles:
        possibles = glob.glob(base_dir + '/*/{}*'.format(species[0]))
    if not possibles:
        possibles = glob.glob(base_dir + '/*/*')
    best = sorted(possibles, key=lambda target: levenshtein_distance(
        species.lower(), construct_species(os.path.basename(target)).lower()))
    return best[0], levenshtein_distance(
        species.lower(), construct_species(os.path.basename(best[0])).lower())


def construct_species(name):
    """Abbreviate capitalized words to one letter"""
    splits = name.split(' ')
    if 'Halobacterium' in splits[0]:
        return name
    new = []
    for i, text in enumerate(splits):
        if text[0] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
            new.append(str(text[0]))
        else:
            new.extend(splits[i:])
            break
    return '_'.join(new)


def read_path_data(path):
    """read the data from the path location"""
    dat = ""
    with open(path) as path_p:
        for next_dat in path_p:
            if next_dat[0] == '>':
                continue
            dat += next_dat.strip()
    return dat


def align(gapped_seq, full_seq):
    """align a gapped sequence to the full sequence"""
    parts = gapped_seq.split('-')
    parts = [p for p in parts if p]  # remove empty strings
    num_parts = len(parts)
    parts_len = sum([len(p) for p in parts])
    full_len = len(full_seq)
    matrix_side_1 = num_parts + 1
    matrix_side_2 = full_len - parts_len + 2
    shape = (matrix_side_1, matrix_side_2)
    costs = np.zeros(shape)
    path = np.zeros_like(costs, dtype=np.int8)
    for i in range(costs.shape[0]):
        costs[i, 0] = np.inf
        path[i, 0] = 1  # take part
    for j in range(costs.shape[1] - 1):
        costs[0, j + 1] = j
        path[0, j + 1] = 2  # take gap
    offset = 0
    for i in range(1, costs.shape[0]):
        part = parts[i - 1]
        size = len(part)
        for j in range(1, costs.shape[1]):
            n = offset + j - 1
            sub_seq = full_seq[n:n+size]
            take_gap = 1 + costs[i, j - 1]
            take_part = len([1 for c1, c2 in zip(sub_seq, part)
                             if c1 != c2])
            take_part += costs[i - 1, j]
            if take_part <= take_gap:
                costs[i, j] = take_part
                path[i, j] = 1  # take part
            else:
                costs[i, j] = take_gap
                path[i, j] = 2  # take gap
        offset += size
    i, j = path.shape
    i -= 1
    j -= 1
    alignment = ""
    while i > 1 or j > 1:
        if path[i, j] == 1:  # take part
            alignment += parts.pop()[::-1]
            i -= 1
        else:  # take gap
            alignment += '-'
            j -= 1
    alignment = alignment[::-1]
    return alignment


if __name__ == "__main__":
    make_no_gaps_to_gaps_file()
