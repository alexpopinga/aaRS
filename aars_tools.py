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
            'amino-acid-paths': predict_amino_path(path, aa),
            'nucleotide-paths': predict_nucleotide_path(path, aa)
        }
        sys.stdout.write('\r{}%'.format((i * 100) // len(no_gaps)))
        sys.stdout.flush()
    with open('no_gaps_to_gaps.json', 'w') as f:
        json.dump(data, f, indent=2)


def predict_amino_path(path, aa):
    """Try to find the amino acid sequence file"""
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


if __name__ == "__main__":
    make_no_gaps_to_gaps_file()
