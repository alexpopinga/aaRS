"""Python tools for aaRS data

a.k.a. "Paul's tools"
"""
import sys
from itertools import product
import json
import xml.etree.ElementTree as ET
import numpy as np

AARS_XML = 'AARS.xml'
BASE_DIR = 'BEAST 2/XMLs/Better Priors (final, actually used XMLs)/'
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


def levenshtein_distance(source, target):
    """Custom Levenshtein distance calculation"""
    distance_matrix = np.zeros((len(source) + 1, len(target) + 1), int)
    distance_matrix[0, :] = np.arange(len(target) + 1)
    distance_matrix[:, 0] = np.arange(len(source) + 1)
    for i, j in product(range(1, len(source)+1), range(1, len(target)+1)):
        substitution_cost = 0 if source[i-1] == target[i-1] else 1
        distance_matrix[i, j] = min(
            distance_matrix[i-1, j] + 1,
            distance_matrix[i, j-1] + 1,
            distance_matrix[i-1, j-1] + substitution_cost
        )
    return int(distance_matrix[len(source), len(target)])


def make_no_gaps_to_gaps_file():
    """Write the no gaps to gaps file"""
    version = 1
    try:
        with open('no_gaps_to_gaps.json') as f:
            data = json.load(f)
        if data['version'] != version:
            data = {'version': version}
    except KeyError:
        data = {'version': version}
    except FileNotFoundError:
        data = {'version': version}
    except json.decoder.JSONDecodeError:
        data = {'version': version}
    no_gaps = final_sequences()
    gaps = all_sequences()
    sys.stdout.write('\t0%')
    sys.stdout.flush()
    for i, no_gaps_key in enumerate(no_gaps.keys()):
        try:
            read_test = {
                'no-gaps-key': data[no_gaps_key]['no-gaps-key'],
                'no-gaps-value': data[no_gaps_key]['no-gaps-value'],
                'gaps-key': data[no_gaps_key]['gaps-key'],
                'gaps-value': data[no_gaps_key]['gaps-value'],
                'distance': data[no_gaps_key]['distance']
            }
        except KeyError:
            gaps_key, distance = find_match(
                no_gaps_key, no_gaps=no_gaps, gaps=gaps)
            data[no_gaps_key] = {
                'no-gaps-key': no_gaps_key,
                'no-gaps-value': no_gaps[no_gaps_key],
                'gaps-key': gaps_key,
                'gaps-value': gaps[gaps_key],
                'distance': distance
            }
            with open('no_gaps_to_gaps.json', 'w') as f:
                json.dump(data, f, indent=2)
        sys.stdout.write('\r\t{}%'.format((i * 100) // len(no_gaps)))
        sys.stdout.flush()


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

    print('\n\nfind_match(seq_name)')
    make_no_gaps_to_gaps_file()
    print()
