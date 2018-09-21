"""Python tools for aaRS data

a.k.a. "Paul's tools"
"""
import xml.etree.ElementTree as ET

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


def is_same_sequence(seq1, seq2):
    """Check if sequences are basically the same

    Ignores gaps, '?', and trailing characters
    """
    if seq1 == "" or seq2 == "":
        return True
    if seq1 == seq2:
        return True
    if seq1[0] == seq2[0]:
        return is_same_sequence(seq1[1:], seq2[1:])
    if seq2[0] == '-' or seq2[0] == '?':
        return is_same_sequence(seq1, seq2[1:])
    if seq2[0] == 'V' and is_same_sequence(seq1, seq2[1:]):
        return True
    return False


def find_matches(seq_name):
    """Find any gapped sequences matching the named sequence

    uses is_same_sequence(seq1, seq2) to find matches
    """
    try:
        seq1 = final_sequences()[seq_name]
    except ValueError:
        print('sequences name', seq_name, 'not found in final sequences file')
        return None
    matches = []
    for seq2_name, seq2 in all_sequences().items():
        if is_same_sequence(seq1, seq2):
            matches.append(seq2_name)
    return matches


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

    print('\n\nis_same_sequence(seq1, seq2)')
    seq1 = 'IIEFSPNIAKPFHAGHLRSTIIGGFLANLIRMNYLGDWGEFSIEKYIDTYARLDVYSGESQLYLTRDVGAAMIYVIASQQDLHAAQFFEILKQMLQHVNFGMSTR'
    seq2 = 'IIEFSPNIAKPFHAGHLRSTIIGGFLANL-IRMNYLGDWG-EFSIEKYIDTYARLD-VYSGESQ-LYLTRDVGAAM-IYVIASQQDLHAAQFFEILKQM-LQHVNF------'
    print('\tseq1 is', seq1)
    print('\tseq2 is', seq2)
    print('\tresult is:', is_same_sequence(seq1, seq2))

    print('\n\nfind_match(seq_name)')
    count_found = 0
    count_not_found = 0
    count_multiple_matches = 0
    for name in final_sequences().keys():
        # print('\tfinding match for:', name, final_sequences()[name])
        found = find_matches(name)
        if found:
            count_found += 1
            if len(found) > 1:
                count_multiple_matches += 1
            # for match in found:
                # print('\t            found:', match, all_sequences()[match], '\n')
        else:
            count_not_found += 1
            # print('\t            found: (no match)\n')
    print('\tfound:', count_found)
    print('\tnot found:', count_not_found)
    print('\tmultiple matches:', count_multiple_matches)
