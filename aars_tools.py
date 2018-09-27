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


def find_match(seq_name, no_gaps=final_sequences(), gaps=all_sequences()):
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
    sys.stdout.write('\t0%')
    sys.stdout.flush()
    data = {}
    for i, no_gaps_key in enumerate(no_gaps.keys()):
        gaps_key, distance = find_match(
            no_gaps_key, no_gaps=no_gaps, gaps=gaps)
        domain, species, aa = parse_no_gaps_keys(no_gaps_key)
        path, path_cost = predict_path(domain, species, aa)
        data[no_gaps_key] = {
            'no-gaps-key': no_gaps_key,
            'domain': domain,
            'species': species,
            'amino acid': aa,
            'gaps-key': gaps_key,
            'no-gaps-value': no_gaps[no_gaps_key],
            'gaps-value': gaps[gaps_key],
            'gaps-sequence-uncertainty': distance,
            'path': path,
            'path-uncertainty': path_cost}
        sys.stdout.write('\r\t{}%'.format((i * 100) // len(no_gaps)))
        sys.stdout.flush()
    with open('no_gaps_to_gaps.json', 'w') as f:
        json.dump(data, f, indent=2)


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


# DEMO
if __name__ == "__main__":

    print('\n\nfinal_sequences() :')
    COUNT = 0
    DATA = final_sequences()
    for k, v in DATA.items():
        COUNT += 1
        if COUNT > 3:
            print('\t...{} more...'.format(len(DATA) - 3))
            break
        print('\t{}: {}'.format(k, v))

    print('\n\nclass_one_sequences() :')
    COUNT = 0
    DATA = class_one_sequences()
    for k, v in DATA.items():
        COUNT += 1
        if COUNT > 3:
            print('\t...{} more...'.format(len(DATA) - 3))
            break
        print('\t{}: {}'.format(k, v))

    print('\n\nclass_two_sequences() :')
    COUNT = 0
    DATA = class_two_sequences()
    for k, v in DATA.items():
        COUNT += 1
        if COUNT > 3:
            print('\t...{} more...'.format(len(DATA) - 3))
            break
        print('\t{}: {}'.format(k, v))

    print('\n\nall_sequences() :')
    COUNT = 0
    DATA = all_sequences()
    for k, v in DATA.items():
        COUNT += 1
        if COUNT > 3:
            print('\t...{} more...'.format(len(DATA) - 3))
            break
        print('\t{}: {}'.format(k, v))

    print('\n\nmake_no_gaps_to_gaps_file() :')
    make_no_gaps_to_gaps_file()
    print()
