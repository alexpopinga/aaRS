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
    data = {}
    for i, no_gaps_key in enumerate(no_gaps.keys()):
        gaps_key, distance = find_match(
            no_gaps_key, no_gaps=no_gaps, gaps=gaps)
        domain, species, aa = parse_no_gaps_keys(no_gaps_key)
        if species == 'R_marinus':
            domain = 'arch'
        if species == 'R_rosetta':
            species = 'S_rosetta'
        path, path_cost = predict_path(domain, species, aa)
        aa_paths = predict_amino_path(path, aa)
        nuc_paths = predict_nucleotide_path(path, aa)
        aa_data = [read_path_data(p) for p in aa_paths]
        aa_align = [align(gaps[gaps_key], d) for d in aa_data]
        aa_visual = [''.join([c1.lower() if c2 == '-' else c1
                              for c1, c2 in zip(dat, ali)])
                     for dat, ali in zip(aa_data, aa_align)]
        aa_misalignments = [sum([1 if c2 not in '-.' and c1 != c2 else 0
                                 for c1, c2 in zip(d, a)])
                            for d, a in zip(aa_data, aa_align)]
        nuc_data = [read_path_data(p) for p in nuc_paths]
        nuc_trans = []
        for nu in nuc_data:
            try:
                aa_idx = [len(a)*3 for a in aa_data].index(len(nu))
            except ValueError:
                nuc_trans.append(translate(nu + ('N' * (3 - len(nu) % 3))))
            else:
                nuc_trans.append(aa_data[aa_idx])
        nuc_trans_align = [align(gaps[gaps_key], d) for d in nuc_trans]
        nuc_misalignments = [sum([1 if c2 not in '-.*' and c1 != c2 else 0
                                  for c1, c2 in zip(d, a)])
                             for d, a in zip(nuc_trans, nuc_trans_align)]
        nuc_align = [''.join(['---' if a[i] in '-.' else n[i*3:(i+1)*3]
                              for i in range(len(a))])
                     for a, n in zip(nuc_trans_align, nuc_data)]
        nuc_visual = [''.join([n[i*3:(i+1)*3].lower() if a[i] in '-.' else n[i*3:(i+1)*3]
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
            'amino-acid-misalignments': aa_misalignments,
            'nucleotide-data': nuc_data,
            'nucleotide-translation': nuc_trans,
            'nucleotide-trans-align': nuc_trans_align,
            'nucleotide-trans-misalignments': nuc_misalignments,
            'nucleotide-aligned-gaps': nuc_align,
            'nucleotide-visual-alignment': nuc_visual
        }
    with open('gap_data.json', 'w') as f:
        json.dump(data, f, indent=2)
    output_split_files()


def output_split_files():
    """Output 2 files with summaries of aligned and unaligned"""
    with open('gap_data.json') as f:
        d = json.load(f)
    perfect = {}  # aligned nucleotide
    good = {}  # aligned amino acids
    bad = {}  # some misaligned amino acids
    really_bad = {}  # very badly aligned
    missing_data = {}  # missing data
    counts = {
        'perfectly aligned to nucleotide': 0,
        'perfectly aligned to amino acids': 0,
        'misaligned to amino acids': 0,
        'very misaligned to amino acids': 0,
        'missing data': 0
    }
    for k, v in d.items():
        mis = v['nucleotide-trans-misalignments']
        mis_alt = v['amino-acid-misalignments']
        if not mis:
            counts['missing data'] += 1
            missing_data[k] = {
                'name': k,
                'regions': v['gaps-value'],
                'best-file': "Missing nucleotide sequence file"
            }
            continue
        if not mis_alt:
            counts['missing data'] += 1
            missing_data[k] = {
                'name': k,
                'regions': v['gaps-value'],
                'best-file': "Missing amino acid sequence file"
            }
            continue
        try:
            idx = mis.index(0)
            counts['perfectly aligned to nucleotide'] += 1
            perfect[k] = {
                'name': k,
                'regions': v['gaps-value'],
                'best-file': v['nucleotide-paths'][idx],
                'aa-num': ''.join(['{0:^3d}'.format(n) for n in range(len(v['nucleotide-translation'][idx]))]),
                'aa-seq': ' ' + '  '.join(list(v['nucleotide-translation'][idx])) + ' ',
                'aa-ali': ' ' + '  '.join(list(v['nucleotide-trans-align'][idx])) + ' ',
                'nuc-al': v['nucleotide-aligned-gaps'][idx]
            }
        except ValueError:
            m = min(mis_alt)
            idx = mis_alt.index(m)
            if m == 0:
                counts['perfectly aligned to amino acids'] += 1
                good[k] = {
                    'name': k,
                    'regions': v['gaps-value'],
                    'best-file': v['amino-acid-paths'][idx],
                    'aa-seq': v['amino-acid-data'][idx],
                    'aa-ali': v['amino-acid-aligned-gaps'][idx]
                }
            elif m < 10:
                counts['misaligned to amino acids'] += 1
                bad[k] = {
                    'name': k,
                    'misalignments': m,
                    'regions': v['gaps-value'],
                    'best-file': v['amino-acid-paths'][idx],
                    'aa-seq': v['amino-acid-data'][idx],
                    'aa-ali': v['amino-acid-aligned-gaps'][idx],
                    'misali': ''.join([' ' if c1 == c2 or c2 in '-.' else '^'
                                       for c1, c2 in zip(v['amino-acid-data'][idx],
                                                         v['amino-acid-aligned-gaps'][idx])])
                }
            else:
                counts['very misaligned to amino acids'] += 1
                really_bad[k] = {
                    'name': k,
                    'misalignments': m,
                    'regions': v['gaps-value'],
                    'best-file': v['amino-acid-paths'][idx],
                    'aa-seq': v['amino-acid-data'][idx],
                    'aa-ali': v['amino-acid-aligned-gaps'][idx],
                    'misali': ''.join([' ' if c1 == c2 or c2 in '-.' else '^'
                                       for c1, c2 in zip(v['amino-acid-data'][idx],
                                                         v['amino-acid-aligned-gaps'][idx])])
                }
    for k, v in counts.items():
        print('{}: {}'.format(k, v))
    with open('gap_data_perfect.txt', 'w') as p:
        json.dump(perfect, p, indent=2)
    with open('gap_data_good.txt', 'w') as g:
        json.dump(good, g, indent=2)
    with open('gap_data_bad.txt', 'w') as b:
        json.dump(bad, b, indent=2)
    with open('gap_data_really_bad.txt', 'w') as rb:
        json.dump(really_bad, rb, indent=2)
    with open('gap_data_missing.txt', 'w') as m:
        json.dump(missing_data, m, indent=2)


def predict_amino_path(path, aa):
    """Try to find the amino acid sequence file"""
    if path[-3:] == "ile" and aa == 'ile':  # M_mobile ends in 'ile' but may not be aa 'ile'
        return predict_amino_path(path, '_ile')
    if "Pyrococcus horikoshii" in path and 'asn' in aa:
        print('Missing amino acid sequence for Pyrococcus horikoshii asn')
        return []
    if "Homo sapiens" in path and 'phe' in aa:
        print('Missing amino acid sequences for Homo sapiens phe')
        return []
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
        raise RuntimeError('Could not find amino acid path for {} ({})'.format(
            os.path.basename(path), aa))
    return paths


def predict_nucleotide_path(path, aa):
    """Try to find the nucleotide sequence file"""
    if path[-3:] == "ile" and aa == 'ile':  # M_modile ends in 'ile' but may not be aa 'ile'
        return predict_nucleotide_path(path, '_ile')
    if 'Methanopyrus kandleri' in path and aa == 'arg':
        print('Corrupt nucleotide sequence for Methanopyrus kandleri arg')
        return []
    if 'Crassostrea gigas' in path and aa == 'ala':
        print('Missing nucleotide sequence for Crassostrea gigas ala')
        return []
    if 'Halobacterium sp.' in path and aa == 'gly':
        print('Corrupt nucleotide sequence for Halobacterium sp. gly')
        return []
    if 'Phycisphaera mikurensis' in path and aa == 'lys':
        print('Missing nucleotide sequence for Phycisphaera mikurensis lys')
        return []
    if "Bos taurus" in path and aa == 'leu':
        print('Missing nucleotide sequence for Bos taurus leu')
        return []
    if "Pyrococcus horikoshii" in path and aa == 'asn':
        print('Missing nucleotide sequence for Pyrococcus horikoshii asn')
        return []
    if "Giardia lamblia" in path and aa == 'asn':
        print('Missing nucleotide sequence for Giardia lamblia asn')
        return []
    if "Homo sapiens" in path and aa == 'phe':
        print('Missing nucleotide sequence for Homo sapiens phe')
        return []
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
    reg = list(gapped_seq.strip('-') + '-')
    seq = list(full_seq)
    gaps = sum([1 if c == '-' else 0 for c in reg])
    matrix_side_1 = len(reg) - gaps + 1
    matrix_side_2 = (len(seq) - (len(reg) - gaps)) + 2
    shape = (matrix_side_1, matrix_side_2)
    costs = np.zeros(shape)
    path = np.zeros_like(costs, dtype=np.int8)
    for i in range(costs.shape[0]):
        costs[i, 0] = np.inf
        path[i, 0] = 0  # invalid
    for j in range(costs.shape[1] - 1):
        costs[0, j + 1] = j
        path[0, j + 1] = 2  # take '-' gap
    i = 0
    for x in range(1, costs.shape[0]):
        c1 = reg[i]
        j = 0
        for y in range(1, costs.shape[1]):
            c2 = seq[j + x - 1]
            gaps_cost = (1 if reg[i + 1] == '-' else
                         10 if path[x, y - 1] != 1 else
                         100
                         )
            take_gap = gaps_cost + costs[x, y - 1]
            take_char = 0 if c1 == c2 else 10000
            take_char += costs[x - 1, y]
            if take_char <= take_gap:
                costs[x, y] = take_char
                path[x, y] = 1  # take part
            else:
                costs[x, y] = take_gap
                path[x, y] = 2 if gaps_cost == 1 else 3  # take '-' or '.'
            j += 1
        i += 1
        while i < len(reg) and reg[i] == '-':
            i += 1

    x, y = path.shape
    x -= 1
    y -= 1
    alignment = []
    reg = [c for c in reg if c != '-']
    while x > 0 or y > 1:
        if path[x, y] == 1:  # take part
            alignment.append(reg.pop())
            x -= 1
        elif path[x, y] == 2:  # take '-' gap
            alignment.append('-')
            y -= 1
        else:  # take gap '.' gap
            alignment.append('.')
            y -= 1
    return ''.join(reversed(alignment))


if __name__ == "__main__":
    make_no_gaps_to_gaps_file()
