"""Python tools for aaRS data

a.k.a. "Paul's tools"
"""
import os.path
import json
import xml.etree.ElementTree as ET
from urllib.parse import urlencode
from urllib.request import urlopen
import re
import tempfile
import zipfile


try:
    from aars_algorithms_fast import levenshtein_distance_c as levenshtein_distance
    from aars_algorithms_fast import align_c as align
except ImportError:
    print('using slow algorithms - this could take a while')
    from aars_algorithms_slow import levenshtein_distance_py as levenshtein_distance
    from aars_algorithms_slow import align_py as align


AARS_XML = 'AARS.xml'
BASE_DIR = 'BEAST 2/XMLs/Better Priors (final, actually used XMLs)/'
A_B_E_ZIP = 'GenBank aa sequences w: no structure/Archaeans_Bacteria_Eukaryotes.zip'
CLASS_I = os.path.join(BASE_DIR, 'ClassI_betterPriors.xml')
CLASS_II = os.path.join(BASE_DIR, 'ClassII_betterPriors.xml')
TMP = tempfile.gettempdir()
OUTPUT_DIR = os.path.join(TMP, 'Archaeans_Bacteria_Eukaryotes/')
LINE_WIDTH = 72

# Result files
RAW_DATA = os.path.join(TMP, 'raw_data.json')


def main():
    """The main code"""
    no_gaps = final_sequences()
    gaps = all_sequences()
    data = {}
    total = len(no_gaps)

    # load existing 'perfect' data
    perfect_keys, existing_data = load_perfect_data()

    for i, no_gaps_key in enumerate(no_gaps.keys()):
        print('({:d}/{:d}) Aligning {}'.format(i+1, total, no_gaps_key))

        # if already perfectly aligned, use that data
        if no_gaps_key in perfect_keys and no_gaps_key in existing_data.keys():
            data[no_gaps_key] = existing_data[no_gaps_key]
            continue

        gaps_key, domain, species, aa, class_type = parse_taxonomy(
            gaps, no_gaps, no_gaps_key)

        ## AMINO ACID DATA ##
        aa_path = os.path.join(
            OUTPUT_DIR, domain, species, '{}.aa'.format(no_gaps_key))
        aa_data = None
        if os.path.exists(aa_path):
            aa_header, aa_id, aa_data, aa_align, aa_misalignments = get_aa_from_file(
                aa_path, gaps, gaps_key)
        if not aa_data:
            aa_header, aa_id, aa_data, aa_align, aa_misalignments = get_aa_from_data_set(
                domain, species, aa, gaps, gaps_key, aa_path)

        ## NUCLEOTIDE DATA ##
        nuc_path = os.path.join(
            OUTPUT_DIR, domain, species, '{}.nuc'.format(no_gaps_key))
        if os.path.exists(nuc_path):
            nuc_header, nuc_id, nuc_data = read_path_file(nuc_path)
        else:
            nuc_header, nuc_id, nuc_data = get_nuc_from_data_set(
                domain, species, aa, nuc_path)

        # create record of data
        data[no_gaps_key] = {
            'aa-gi': aa_id,
            'nuc-gi': nuc_id,
            'nuc_header': nuc_header,
            'aa-header': aa_header,
            'no-gaps-key': no_gaps_key,
            'class': class_type,
            'domain': domain,
            'species': species,
            'amino-acid': aa,
            'gaps-key': gaps_key,
            'no-gaps-value': no_gaps[no_gaps_key],
            'gaps-value': gaps[gaps_key],
            'amino-acid-path': aa_path,
            'nucleotide-path': nuc_path,
            'amino-acid-data': aa_data,
            'amino-acid-aligned-gaps': aa_align,
            'amino-acid-misalignments': aa_misalignments,
            'nucleotide-data': nuc_data
        }
    with open(RAW_DATA, 'w') as file_p:
        json.dump(data, file_p, indent=2)
    output_split_files()


def final_sequences():
    """dictionary of final sequences (no gaps)"""
    tree = ET.parse(AARS_XML)
    root = tree.getroot()
    data = root[0]
    return {
        s.attrib['taxon']: s.attrib['value']
        for s in data if s.tag == 'sequence'
        and '1f7u' not in s.attrib['taxon']
        and '1iq0' not in s.attrib['taxon']
        and '2zue' not in s.attrib['taxon']
        and '3fnr' not in s.attrib['taxon']
        and '1li5' not in s.attrib['taxon']
        and '1i6m' not in s.attrib['taxon']
        and '1r6u' not in s.attrib['taxon']
        and '2ip1' not in s.attrib['taxon']
        and '2dlc' not in s.attrib['taxon']
        and '1x54a' not in s.attrib['taxon']
        and '3m4p' not in s.attrib['taxon']
        and '3i7f' not in s.attrib['taxon']
        and '2zt5' not in s.attrib['taxon']
        and '3hri' not in s.attrib['taxon']
        and '3lc0' not in s.attrib['taxon']
        and '3bju' not in s.attrib['taxon']
        and '3l4g' not in s.attrib['taxon']
        and '1wleb' not in s.attrib['taxon']
        and 'd1e1oa' not in s.attrib['taxon']
        and 'd1sera' not in s.attrib['taxon']
        and '3hxv' not in s.attrib['taxon']
        and '2zp1' not in s.attrib['taxon']
        and '1wq4' not in s.attrib['taxon']
        and '2o5r' not in s.attrib['taxon']
        and '1qtq' not in s.attrib['taxon']
        and '1nzj' not in s.attrib['taxon']
        and '2cfo' not in s.attrib['taxon']
        and '1rqg' not in s.attrib['taxon']
        and '2cya' not in s.attrib['taxon']
        and 'd1asza' not in s.attrib['taxon']
        and '1wydb' not in s.attrib['taxon']
        and 'd1b8aa' not in s.attrib['taxon']
        and '3g1z' not in s.attrib['taxon']
        and '2rhq' not in s.attrib['taxon']
        and '2i4la' not in s.attrib['taxon']
        and '2j3mb' not in s.attrib['taxon']
        and '2cjab' not in s.attrib['taxon']
        and '3a32' not in s.attrib['taxon']
        and '3c8z' not in s.attrib['taxon']
        and '1htt' not in s.attrib['taxon']
        and '1wu7a' not in s.attrib['taxon']
        and 'd1adjc' not in s.attrib['taxon']
        and 'd1qe0a' not in s.attrib['taxon']
        and '2hz7' not in s.attrib['taxon']
        and '1j09' not in s.attrib['taxon']
        and '1ile' not in s.attrib['taxon']
        and '1qu2' not in s.attrib['taxon']
        and '1wz2' not in s.attrib['taxon']
        and '2v0c' not in s.attrib['taxon']
        and '1pfv' not in s.attrib['taxon']
        and '1woy' not in s.attrib['taxon']
        and '2csx' not in s.attrib['taxon']
        and '2x1l' not in s.attrib['taxon']
        and '2el7' not in s.attrib['taxon']
        and '2g36' not in s.attrib['taxon']
        and '1h3f' not in s.attrib['taxon']
        and '1jil' not in s.attrib['taxon']
        and '2cyb' not in s.attrib['taxon']
        and '2cyc' not in s.attrib['taxon']
        and '1yfsa' not in s.attrib['taxon']
        and '2ztg' not in s.attrib['taxon']
        and '2zze' not in s.attrib['taxon']
        and 'd1c0aa' not in s.attrib['taxon']
        and 'd1efwa' not in s.attrib['taxon']
        and 'd1atib' not in s.attrib['taxon']
        and '1htt' not in s.attrib['taxon']
        and '1wu7a' not in s.attrib['taxon']
        and 'd1adjc' not in s.attrib['taxon']
        and 'd1qe0a' not in s.attrib['taxon']
        and '3e9h' not in s.attrib['taxon']
        and '2alya' not in s.attrib['taxon']
        and 'd1hc7a' not in s.attrib['taxon']
        and 'd1nj8a' not in s.attrib['taxon']
        and 'd1evka' not in s.attrib['taxon']
        and 'd1nyra' not in s.attrib['taxon']
        and '1gax' not in s.attrib['taxon']
        and '3ial' not in s.attrib['taxon']
    }


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
    return best_key


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


def parse_taxonomy(gaps, no_gaps, no_gaps_key):
    '''parse taxonomy details from region key'''
    gaps_key = find_match(no_gaps_key, no_gaps=no_gaps, gaps=gaps)
    domain, species, aa = parse_no_gaps_keys(no_gaps_key)
    class_type = "class_I" if gaps_key in class_one_sequences().keys() else "class_II"

    # TODO Is R. marinus in the wrong domain?
    # For now, I'm manually changing the domain to archaea.
    if species == 'R_marinus':
        domain = 'Archaea'
    ##################################################################

    # TODO is R. rosetta correct?
    # I'm assuming this is actually S. rosetta.
    if species == 'R_rosetta':
        species = 'S_rosetta'
    ##################################################################

    return gaps_key, domain, species, aa, class_type


def load_perfect_data():
    """load the data from any previous run that was already perfect"""
    try:
        with open(os.path.join(TMP, 'gap_data_perfect.txt')) as keys_file:
            perfect_keys = json.load(keys_file).keys()
    except FileNotFoundError:
        perfect_keys = {}
    try:
        with open(RAW_DATA) as data_file:
            existing_data = json.load(data_file)
    except FileNotFoundError:
        existing_data = {}
    return perfect_keys, existing_data


def get_nuc_from_data_set(domain, species, aa, cache_path):
    """read nuc data from Alex's source data set (the zip file)"""
    path, _ = predict_path(domain, species, aa)
    nuc_paths = predict_nucleotide_path(path, aa)
    nuc_paths = [p for p in nuc_paths if 
        "Aaeolicus_arg_nuc(2)" not in p
        and "Mdomestica_val_nuc copy" not in p]

    # TODO What is this file?
    # Archaeans_Bacteria_Eukaryotes/Bacteria/Chroococcidiopsis thermalis/amino acid sequences/Cthermalis_asn_asp_nuc
    # Need to fix and remove exception code.
    nuc_paths = [
        p for p in nuc_paths if 'Cthermalis_asn_asp_nuc' not in p]
    ###############################################################################################################

    if len(nuc_paths) > 1:
        raise RuntimeError(
            'Too many paths found!\n{}'.format(nuc_paths))
    try:
        nuc_path = nuc_paths[0]
    except IndexError:
        nuc_path = None
    with zipfile.ZipFile(A_B_E_ZIP) as z_file:
        nuc_header, nuc_id, nuc_data = read_zip_file_data(nuc_path, z_file)
    # cache this data
    directory = os.path.dirname(cache_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(cache_path, 'w') as file_p:
        file_p.write('>gi|' + (nuc_id if nuc_id else 'none') + '\n' + nuc_data + '\n')
    return nuc_header, nuc_id, nuc_data


def get_aa_from_data_set(domain, species, aa, gaps, gaps_key, cache_path):
    """read aa data from Alex's source data set (the zip file)"""
    path, _ = predict_path(domain, species, aa)
    aa_paths = [p for p in predict_amino_path(path, aa)
                if "CUT" not in p and "TEST" not in p]

    # TODO What is this file?
    # Archaeans_Bacteria_Eukaryotes/Bacteria/Chroococcidiopsis thermalis/amino acid sequences/Cthermalis_asn_asp_aa
    # Need to fix and remove exception code.
    aa_paths = [
        p for p in aa_paths if 'Cthermalis_asn_asp_aa' not in p]
    ###############################################################################################################

    if len(aa_paths) > 1:
        raise RuntimeError(
            'Too many paths found!\n{}'.format(aa_paths))
    try:
        aa_path = aa_paths[0]
    except IndexError:
        print('{} {} {} amino acid data missing from source data set'.format(
            domain, species, aa
        ))
        aa_path = None
    with zipfile.ZipFile(A_B_E_ZIP) as z_file:
        aa_header, aa_id, aa_data = read_zip_file_data(aa_path, z_file)
        if aa_data:
            aa_align = align(gaps[gaps_key], aa_data)
            aa_misalignments = sum([1 if c2 not in '-.*?' and c1 not in '*?' and c1 != c2 else 0
                                    for c1, c2 in zip(aa_data, aa_align)])
        else:
            aa_align = None
            aa_misalignments = None
    if not aa_data:
        raise RuntimeError('No data for {} in source data set'.format(cache_path))
    # cache this data
    directory = os.path.dirname(cache_path)
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(cache_path, 'w') as file_p:
        file_p.write('>gi|' + (aa_id if aa_id else 'none') + '\n' + aa_data + '\n')
    return aa_header, aa_id, aa_data, aa_align, aa_misalignments


def get_aa_from_file(aa_path, gaps, gaps_key):
    """read from a file to get the data"""
    aa_header, aa_id, aa_data = read_path_file(aa_path)
    if aa_data:
        aa_align = align(gaps[gaps_key], aa_data)
        aa_misalignments = sum([1 if c2 not in '-.*?' and c1 not in '*?' and c1 != c2 else 0
                                for c1, c2 in zip(aa_data, aa_align)])
    else:
        aa_align = None
        aa_misalignments = None
    return aa_header, aa_id, aa_data, aa_align, aa_misalignments


def output_split_files():
    """Output files with summaries of aligned and unaligned"""
    with open(RAW_DATA) as file_p:
        d = json.load(file_p)
    perfect = {}  # aligned nucleotide
    genbank = {}  # aligned to GenBank data
    good = {}  # aligned amino acids
    bad = {}  # some misaligned amino acids
    missing_data = {}  # missing data
    counts = {
        'perfectly aligned to our nucleotide data': 0,
        'perfectly aligned to GenBank nucleotide data': 0,
        'aligned to amino acids only': 0,
        'misaligned to amino acids': 0,
        'missing data': 0
    }
    total = len(d)
    for current, (k, v) in enumerate(d.items()):
        print('({:d}/{:d}) Checking {}'.format(current+1, total, k))
        if not v['nucleotide-path']:
            if not v['aa-gi'] or not check_genbank(v['aa-gi'])[0]:
                counts['missing data'] += 1
                missing_data[k] = {
                    'name': k,
                    'comment': "Missing nucleotide sequence file"
                }
                continue
        if not v['amino-acid-path']:
            counts['missing data'] += 1
            missing_data[k] = {
                'name': k,
                'comment': "Missing amino acid sequence file"
            }
            continue
        misalignment = v['amino-acid-misalignments']
        estimated_nuc_len = len(v['amino-acid-data']) * 3
        actual_nuc_len = len(v['nucleotide-data'])
        if misalignment == 0 and (estimated_nuc_len == actual_nuc_len
                                  # From Peter Wills - 10 Oct 2018
                                  # "Many of these aaRS sequences are of the
                                  # 'catalytic domain' and some have other
                                  # domains attached to them, so you may have
                                  # the first codon of the following domain."
                                  #
                                  # Given this, if the nucleotide sequences has
                                  # one extra codon, we ignore it.
                                  or actual_nuc_len - estimated_nuc_len <= 3):

            counts['perfectly aligned to our nucleotide data'] += 1
            nuc_align = ''.join(['---' if v['amino-acid-aligned-gaps'][i] in '-.'
                                 else v['nucleotide-data'][i*3:(i+1)*3]
                                 for i in range(len(v['amino-acid-aligned-gaps']))])
            perfect[k] = {
                'name': k,
                'class': v['class'],
                'regions': v['gaps-value'],
                'amino-acid-file': v['amino-acid-path'],
                'nucleotide-file': v['nucleotide-path'],
                'aa-num': ''.join(['{0:^3d}'.format(n) for n in range(1, len(v['amino-acid-data']) + 1)]),
                'aa-seq': ' ' + '  '.join(list(v['amino-acid-data'])) + ' ',
                'aa-ali': ' ' + '  '.join(list(v['amino-acid-aligned-gaps'])) + ' ',
                'nuc-al': nuc_align
            }
            continue
        # see if GenBank can do better
        # currently, we only check GenBank for the GI parsed from the amino acid
        # file and ignore the GI from the nucleotide file
        check, gb_aa, gb_nuc = check_genbank(v['aa-gi'])
        if check:
            aa_align = align(v['gaps-value'], gb_aa)
            aa_misalignments = sum([1 if c2 not in '-.*?' and c1 not in '*?' and c1 != c2 else 0
                                    for c1, c2 in zip(gb_aa, aa_align)])
            if aa_misalignments == 0:
                counts['perfectly aligned to GenBank nucleotide data'] += 1
                nuc_align = ''.join(['---' if v['amino-acid-aligned-gaps'][i] in '-.'
                                     else gb_nuc[i*3:(i+1)*3]
                                     for i in range(len(v['amino-acid-aligned-gaps']))])
                genbank[k] = {
                    'name': k,
                    'class': v['class'],
                    'regions': v['gaps-value'],
                    'amino-acid-file': 'downloaded from GenBank: gi|{}'.format(v['aa-gi']),
                    'nucleotide-file': 'downloaded from GenBank: gi|{}'.format(v['aa-gi']),
                    'aa-num': ''.join(['{0:^3d}'.format(n) for n in range(1, len(gb_aa) + 1)]),
                    'aa-seq': ' ' + '  '.join(list(gb_aa)) + ' ',
                    'aa-ali': ' ' + '  '.join(list(aa_align)) + ' ',
                    'nuc-al': nuc_align
                }

                # replace data in cached copy of data set
                aa_cache_path = os.path.join(
                    OUTPUT_DIR, v['domain'], v['species'], '{}.aa'.format(k)
                )
                nuc_cache_path = os.path.join(
                    OUTPUT_DIR, v['domain'], v['species'], '{}.nuc'.format(k)
                )
                directory = os.path.dirname(aa_cache_path)
                if not os.path.exists(directory):
                    os.makedirs(directory)
                with open(aa_cache_path, 'w') as file_p:
                    file_p.write('>gi|' + v['aa-gi'] + '|hello\n' + gb_aa + '\n')
                with open(nuc_cache_path, 'w') as file_p:
                    file_p.write('>gi|' + v['aa-gi'] + '|hello\n' + gb_nuc + '\n')
                continue
        if misalignment == 0:
            if estimated_nuc_len < actual_nuc_len:
                too_long = actual_nuc_len - estimated_nuc_len
                comment = 'nucleotide data is {} characters too long'.format(
                    too_long)
            else:
                too_short = estimated_nuc_len - actual_nuc_len
                comment = 'nucleotide data is {} characters too short'.format(
                    too_short)
            counts['aligned to amino acids only'] += 1
            good[k] = {
                'name': k,
                'class': v['class'],
                'aa-genbank-id': v['aa-gi'],
                'nuc-genbank-id': v['nuc-gi'],
                'regions': v['gaps-value'],
                'amino-acid-file': v['amino-acid-path'],
                'nucleotide-file': v['nucleotide-path'],
                'aa-seq': v['amino-acid-data'],
                'aa-ali': v['amino-acid-aligned-gaps'],
                'comment': comment
            }
            continue
        counts['misaligned to amino acids'] += 1
        bad[k] = {
            'name': k,
            'class': v['class'],
            'misalignments': misalignment,
            'regions': v['gaps-value'],
            'amino-acid-file': v['amino-acid-path'],
            'nucleotide-file': v['nucleotide-path'],
            'aa-seq': v['amino-acid-data'],
            'aa-ali': v['amino-acid-aligned-gaps'],
            'misali': ''.join([' ' if c1 == c2 or c1 in '*?' or c2 in '-.*?' else '^'
                               for c1, c2 in zip(v['amino-acid-data'],
                                                 v['amino-acid-aligned-gaps'])])
        }
    # write JSON output data
    with open(os.path.join(TMP, 'gap_data_perfect.txt'), 'w') as p:
        json.dump(perfect, p, indent=2)
    with open(os.path.join(TMP, 'gap_data_genbank.txt'), 'w') as gb:
        json.dump(genbank, gb, indent=2)
    with open(os.path.join(TMP, 'gap_data_good.txt'), 'w') as g:
        json.dump(good, g, indent=2)
    with open(os.path.join(TMP, 'gap_data_bad.txt'), 'w') as b:
        json.dump(bad, b, indent=2)
    with open(os.path.join(TMP, 'gap_data_missing.txt'), 'w') as m:
        json.dump(missing_data, m, indent=2)

    # write perfect HTML output data
    with open('gap_data_perfect.html', 'w') as file_p:
        file_p.write('<h1>Perfectly aligned sequences</h1>\n')
        # HTML header
        file_p.write('<!DOCTYPE html>\n')
        file_p.write('<html>\n')
        file_p.write('<head>\n')
        file_p.write('<title>Perfectly aligned sequences</title>\n')
        file_p.write('</head>\n\n')

        # HTML body
        file_p.write('<body>\n\n')

        for _, aa in perfect.items():
            file_p.write('<h2>{}</h1>\n'.format(aa['name']))
            file_p.write('<p>Class {}</p>\n'.format(
                'I' if aa['class'] == 'class_I' else 'II'
            ))
            file_p.write('<p>{}</p>\n'.format(
                aa['amino-acid-file']
            ))
            file_p.write('<p>{}</p>\n'.format(
                aa['nucleotide-file']
            ))

            file_p.write('<h3>alignment summary</h3>\n')
            file_p.write('<p style="font-family:monospace">\n')
            for i in range(0, len(aa['nuc-al']), LINE_WIDTH):
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-num'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-num'][j])
                    except IndexError:
                        break
                file_p.write('<br/>\n')
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-seq'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-seq'][j])
                    except IndexError:
                        break
                file_p.write('<br/>\n')
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-ali'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-ali'][j])
                    except IndexError:
                        break
                file_p.write('<br/>\n')
                for j in range(i, i + LINE_WIDTH):
                    try:
                        file_p.write(aa['nuc-al'][j])
                    except IndexError:
                        break
                file_p.write('<br/><br/>\n')
            file_p.write('</p>\n\n')

        file_p.write('</body>\n')
        file_p.write('</html>\n')

    # write GenBank aligned HTML output data
    with open('gap_data_genbank.html', 'w') as file_p:
        file_p.write('<h1>GenBank aligned sequences</h1>\n')
        # HTML header
        file_p.write('<!DOCTYPE html>\n')
        file_p.write('<html>\n')
        file_p.write('<head>\n')
        file_p.write('<title>GenBank aligned sequences</title>\n')
        file_p.write('</head>\n\n')

        # HTML body
        file_p.write('<body>\n\n')

        for _, aa in genbank.items():
            file_p.write('<h2>{}</h1>\n'.format(aa['name']))
            file_p.write('<p>Class {}</p>\n'.format(
                'I' if aa['class'] == 'class_I' else 'II'
            ))
            file_p.write('<p>{}</p>\n'.format(
                aa['amino-acid-file']
            ))
            file_p.write('<p>{}</p>\n'.format(
                aa['nucleotide-file']
            ))

            file_p.write('<h3>alignment summary</h3>\n')
            file_p.write('<p style="font-family:monospace">\n')
            for i in range(0, len(aa['nuc-al']), LINE_WIDTH):
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-num'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-num'][j])
                    except IndexError:
                        break
                file_p.write('<br/>\n')
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-seq'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-seq'][j])
                    except IndexError:
                        break
                file_p.write('<br/>\n')
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-ali'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-ali'][j])
                    except IndexError:
                        break
                file_p.write('<br/>\n')
                for j in range(i, i + LINE_WIDTH):
                    try:
                        file_p.write(aa['nuc-al'][j])
                    except IndexError:
                        break
                file_p.write('<br/><br/>\n')
            file_p.write('</p>\n\n')

        file_p.write('</body>\n')
        file_p.write('</html>\n')

    # write amino acid aligned HTML output data
    with open('gap_data_good.html', 'w') as file_p:
        file_p.write('<h1>Not aligned to nucleotide data</h1>\n')
        # HTML header
        file_p.write('<!DOCTYPE html>\n')
        file_p.write('<html>\n')
        file_p.write('<head>\n')
        file_p.write('<title>Not aligned to nucleotide data</title>\n')
        file_p.write('</head>\n\n')

        # HTML body
        file_p.write('<body>\n\n')

        for _, aa in good.items():
            file_p.write('<h2>{}</h1>\n'.format(aa['name']))
            file_p.write('<p>Class {}</p>\n'.format(
                'I' if aa['class'] == 'class_I' else 'II'
            ))
            file_p.write('<p>{}</p>\n'.format(
                aa['amino-acid-file']
            ))
            file_p.write('<p>{}</p>\n'.format(
                aa['nucleotide-file']
            ))
            file_p.write('<p>{}</p>\n'.format(
                aa['comment']
            ))

            file_p.write('<h3>alignment summary</h3>\n')
            file_p.write('<p style="font-family:monospace">\n')
            for i in range(0, len(aa['aa-seq']), LINE_WIDTH):
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-seq'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-seq'][j])
                    except IndexError:
                        break
                file_p.write('<br/>\n')
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-ali'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-ali'][j])
                    except IndexError:
                        break
                file_p.write('<br/><br/>\n')
            file_p.write('</p>\n\n')

        file_p.write('</body>\n')
        file_p.write('</html>\n')

    # write bad HTML output data
    with open('gap_data_bad.html', 'w') as file_p:
        file_p.write('<h1>Misaligned to amino acid sequences</h1>\n')
        # HTML header
        file_p.write('<!DOCTYPE html>\n')
        file_p.write('<html>\n')
        file_p.write('<head>\n')
        file_p.write('<title>Misaligned to amino acid sequences</title>\n')
        file_p.write('</head>\n\n')

        # HTML body
        file_p.write('<body>\n\n')

        for _, aa in bad.items():
            file_p.write('<h2>{}</h1>\n'.format(aa['name']))
            file_p.write('<p>Class {}</p>\n'.format(
                'I' if aa['class'] == 'class_I' else 'II'
            ))
            file_p.write('<p>{}</p>\n'.format(
                aa['amino-acid-file']
            ))
            file_p.write('<p>{} misalignment{}</p>\n'.format(
                aa['misalignments'],
                '' if aa['misalignments'] == 1 else 's'
            ))

            file_p.write('<h3>alignment summary</h3>\n')
            file_p.write('<p style="font-family:monospace">\n')
            for i in range(0, len(aa['aa-seq']), LINE_WIDTH):
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-seq'][j] == ' ':
                            file_p.write('&nbsp;')
                        else:
                            file_p.write(aa['aa-seq'][j])
                    except IndexError:
                        break
                file_p.write('<br/>\n')
                for j in range(i, i + LINE_WIDTH):
                    try:
                        if aa['aa-ali'][j] == ' ':
                            c = '&nbsp;'
                        else:
                            c = aa['aa-ali'][j]
                        if aa['misali'][j] == '^':
                            file_p.write(bad_style(c))
                        else:
                            file_p.write(c)
                    except IndexError:
                        break
                file_p.write('<br/><br/>\n')
            file_p.write('</p>\n\n')

        file_p.write('</body>\n')
        file_p.write('</html>\n')

    # write missing HTML output data
    with open('gap_data_missing.html', 'w') as file_p:
        file_p.write('<h1>Missing data</h1>\n')
        # HTML header
        file_p.write('<!DOCTYPE html>\n')
        file_p.write('<html>\n')
        file_p.write('<head>\n')
        file_p.write('<title>Missing data</title>\n')
        file_p.write('</head>\n\n')

        # HTML body
        file_p.write('<body>\n\n')

        for _, aa in missing_data.items():
            file_p.write('<h2>{}</h1>\n'.format(aa['name']))
            file_p.write('<p>{}</p>\n'.format(
                aa['comment']
            ))

        file_p.write('</body>\n')
        file_p.write('</html>\n')

    # write CSV class I output data
    with open(os.path.join(TMP, 'gap_data_alignment_class_I.csv'), 'w') as csv_file:
        csv_data = {k: v for k, v in perfect.items() if v['class'] == 'class_I'}
        write_csv(d, csv_file, csv_data)

    # write CSV class II output data
    with open(os.path.join(TMP, 'gap_data_alignment_class_II.csv'), 'w') as csv_file:
        csv_data = {k: v for k, v in perfect.items() if v['class'] == 'class_II'}
        write_csv(d, csv_file, csv_data)

    # print summary data
    for k, v in counts.items():
        print('{}: {}'.format(k, v))


def write_csv(d, csv_file, csv_data):
    '''write CSV data'''
    header = 'Identifier,Amino Acid'
    counter = 1
    total = len(csv_data)
    indices = [0 for _ in csv_data]
    strings = ['{},{}'.format(
        k, d[k]['amino-acid']) for k in csv_data]
    gap = True

    while any(i < len(d[k]['amino-acid-aligned-gaps']) for i, k in zip(indices, csv_data)):
        header += ',{}'.format(counter)
        counter += 1

        # first, check if in a non-gap area with a removed section (i.e. a '.' character)
        dot = not gap and any(
            i < len(d[k]['amino-acid-aligned-gaps']) and d[k]['amino-acid-aligned-gaps'][i] == '.' for i, k in zip(indices, csv_data))

        for i, k in enumerate(csv_data):
            v = d[k]

            # check if at end of alignment
            if indices[i] >= len(v['amino-acid-aligned-gaps']):
                continue  # no need to write anything

            # get character
            c = v['amino-acid-aligned-gaps'][indices[i]]

            strings[i] += ','
            if gap and c == '-':  # only '-' alignments may advance
                strings[i] += '-'
                indices[i] += 1
                continue
            elif dot and c == '.':  # only '.' alignments may advance
                strings[i] += '.'
                indices[i] += 1
                continue
            elif not gap and c not in '.-':
                strings[i] += str(indices[i] + 1)
                indices[i] += 1
                continue

        # if all indices are not '-', we change to non-gap
        if all(i >= len(d[k]['amino-acid-aligned-gaps']) or d[k]['amino-acid-aligned-gaps'][i] != '-' for i, k in zip(indices, csv_data)):
            gap = False
        # if all indices are '-', we change to gap
        if all(i >= len(d[k]['amino-acid-aligned-gaps']) or d[k]['amino-acid-aligned-gaps'][i] == '-' for i, k in zip(indices, csv_data)):
            gap = True
        # otherwise, gap status stays the same

    csv_file.write(header + '\n')
    csv_file.write('\n'.join(strings) + '\n')


def predict_amino_path(path, aa):
    """Try to find the amino acid sequence file"""
    file1 = re.compile(r'^{}/amino[^/]*/[^/]*{}_[^/]*$'.format(path, aa))
    file2 = re.compile(r'^{}/*/[^/]*{}_aa[^/]*$'.format(path, aa))

    with zipfile.ZipFile(A_B_E_ZIP) as z_file:
        files = z_file.namelist()

    if path[-3:] == "ile" and aa == 'ile':  # M_mobile ends in 'ile' but may not be aaRS 'ile'
        return predict_amino_path(path, '_ile')
    paths = list(filter(file1.match, files))
    if not paths:
        paths = list(filter(file2.match, files))
    if not paths:
        if aa == 'glu':  # glu is sometimes listed as gluPro
            return predict_amino_path(path, 'gluPro')
        if aa == 'leu':  # leu is sometimes listed as leuALPHA
            return predict_amino_path(path, 'leuALPHA')
        if aa == 'met':  # met is sometimes listed as metALPHA
            return predict_amino_path(path, 'metALPHA')
        if aa == 'phe':  # phe is sometimes listed as pheALPHA
            return predict_amino_path(path, 'pheALPHA')
        if aa == 'tyr':  # tyr is sometimes listed as Tyr
            return predict_amino_path(path, 'Tyr')
        if aa == 'val':  # val is sometimes listed as Val or val1
            try:
                paths = predict_amino_path(path, 'Val')
            except RuntimeError:
                paths = predict_amino_path(path, 'val1')
    return paths


def predict_nucleotide_path(path, aa):
    """Try to find the nucleotide sequence file"""
    file1 = re.compile(r'^{}/nuc[^/]*/[^/]*{}_[^/]*$'.format(path, aa))
    file2 = re.compile(r'^{}/*/[^/]*{}_nuc[^/]*$'.format(path, aa))

    with zipfile.ZipFile(A_B_E_ZIP) as z_file:
        files = z_file.namelist()

    if path[-3:] == "ile" and aa == 'ile':  # M_modile ends in 'ile' but may not be aaRS 'ile'
        return predict_nucleotide_path(path, '_ile')
    paths = list(filter(file1.match, files))
    if not paths:
        paths = list(filter(file2.match, files))
    if not paths:
        if aa == 'glu':  # glu is sometimes listed as gluPro
            return predict_nucleotide_path(path, 'gluPro')
        if aa == 'leu':  # leu is sometimes listed as leuALPHA
            return predict_nucleotide_path(path, 'leuALPHA')
        if aa == 'met':  # met is sometimes listed as metALPHA
            return predict_nucleotide_path(path, 'metALPHA')
        if aa == 'phe':  # phe is sometimes listed as pheALPHA
            return predict_nucleotide_path(path, 'pheALPHA')
    return paths


def predict_path(domain, species, aa):
    """Try to find the path to the original sequences"""
    path1 = re.compile(r'^{}/{}[^/]*/$'.format(domain, species[0]))
    path2 = re.compile(r'^{}/[^/]*/$'.format(domain))
    path3 = re.compile(
        r'^(Archaea|Bacteria|Eukaryote)/{}[^/]*/$'.format(species[0]))
    path4 = re.compile(r'^(Archaea|Bacteria|Eukaryote)/[^/]*/$')

    with zipfile.ZipFile(A_B_E_ZIP) as z_file:
        files = z_file.namelist()

    possibles = list(filter(path1.match, files))
    if not possibles:
        possibles = list(filter(path2.match, files))
    if not possibles:
        possibles = list(filter(path3.match, files))
    if not possibles:
        possibles = list(filter(path4.match, files))
    possibles = [p[:-1] for p in possibles]

    best = sorted(possibles, key=lambda target: levenshtein_distance(
        species.lower(), construct_species(os.path.basename(target)).lower()))
    return best[0], levenshtein_distance(
        species.lower(), construct_species(os.path.basename(best[0])).lower())


def construct_species(name):
    """Abbreviate capitalized words to one letter"""
    if name == '':
        return '_'
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


def read_zip_file_data(path, z_file):
    """read the data from the zip file"""
    if path is None:
        return None, None, None
    dat = ''
    header = None
    gi = None
    with z_file.open(path) as path_p:
        for next_dat_b in path_p.readlines():
            next_dat = next_dat_b.decode()
            if next_dat[0] == '>':
                header = next_dat.strip()
                if next_dat[0:4] == '>gi|':
                    try:
                        gi = next_dat[4:].split('|')[0].split()[0]
                    except IndexError:
                        gi = None
                else:
                    gi = None
                continue
            dat += next_dat.strip()
    return header, gi, dat


def read_path_file(path):
    """read the data from the path file"""
    if path is None:
        return None, None, None
    dat = ''
    header = None
    gi = None
    with open(path) as path_p:
        for next_dat in path_p.readlines():
            if next_dat[0] == '>':
                header = next_dat.strip()
                if next_dat[0:4] == '>gi|':
                    try:
                        gi = next_dat[4:].split('|')[0].split()[0]
                    except IndexError:
                        gi = None
                else:
                    gi = None
                continue
            dat += next_dat.strip()
    return header, gi, dat


def check_genbank(gi):
    """check if genbank data could be used"""
    if not gi:
        return False, None, None
    # check for cached nuc data
    try:
        with open(os.path.join(TMP, '{}.nuc'.format(gi))) as nuc_file:
            nuc_data = nuc_file.read()
            nuc_len = len(nuc_data)
    except FileNotFoundError:
        # no cached nuc data -- re-download
        nuc_data = get_genbank_nuc(gi)
        nuc_len = len(nuc_data)
        with open(os.path.join(TMP, '{}.nuc'.format(gi)), 'w') as nuc_file:
            nuc_file.write(nuc_data)

    # check for cached aa data
    try:
        with open(os.path.join(TMP, '{}.aa'.format(gi))) as aa_file:
            aa_data = aa_file.read()
            aa_len = 3 * len(aa_data)
    except FileNotFoundError:
        # no cached aa data -- re-download
        aa_data = get_genbank_aa(gi)
        aa_len = 3 * len(aa_data)
        with open(os.path.join(TMP, '{}.aa'.format(gi)), 'w') as aa_file:
            aa_file.write(aa_data)

    if nuc_len > 6 and (aa_len == nuc_len or aa_len == nuc_len - 3):
        return True, aa_data, nuc_data
    return False, aa_data, nuc_data


def get_genbank_aa(gi):
    """get aa data from GenBank"""
    return get_genbank(gi, 'aa')


def get_genbank_nuc(gi):
    """get aa data from GenBank"""
    return get_genbank(gi, 'na')


def get_genbank(gi, type_):
    """get data from GenBank"""
    try:
        gi_, subseq = gi.split(':')
        start, stop = subseq.split('-')
        rettype = 'fasta_cds_' + type_
        params = urlencode({
            'db': 'nuccore',
            'id': gi_,
            'seq_start': start,
            'seq_stop': stop,
            'rettype': rettype,
            'retmode': 'text'
        })
    except ValueError:
        rettype = 'fasta_cds_' + type_
        params = urlencode({
            'db': 'nuccore',
            'id': gi,
            'rettype': rettype,
            'retmode': 'text'
        })
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + params
    dat = ''
    print(url)
    with urlopen(url) as web:
        for line in web.read().decode('UTF-8').splitlines():
            if line and line[0] == '>':
                continue
            for c in line:
                if c == '\n':
                    continue
                if c not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ-*':
                    return None
                dat += c
    return dat


def pretty_print_alignment(seq, align, width=72):
    """pretty print the alignment"""
    i = 0
    while i < len(seq):
        print('{}\n{}\n\n'.format(seq[i:i+width], align[i:i+width]))
        i += width


def bad_style(c):
    """highlight a span of text"""
    return '<span style="background-color:yellow">'+c+'</span>'


def load_missing_data():
    """manually enter data into temp directory"""
    # TODO Add proper aa data for leu_arch_T_volcanium into source data set
    directory = os.path.join(OUTPUT_DIR, 'Archaea', 'T_volcanium')
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename = os.path.join(directory, 'leu_arch_T_volcanium.aa')
    with open(filename, 'w') as f:
        f.write(LEU_ARCH_T_VOLCANIUM_AA)
    #######################################################################

    # TODO Add proper data for trp_arch_T_volcanium into source data set
    directory = os.path.join(OUTPUT_DIR, 'Archaea', 'T_volcanium')
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename_aa = os.path.join(directory, 'trp_arch_T_volcanium.aa')
    with open(filename_aa, 'w') as f:
        f.write(TRP_ARCH_T_VOLCANIUM_AA)
    filename_nuc = os.path.join(directory, 'trp_arch_T_volcanium.nuc')    
    with open(filename_nuc, 'w') as f:
        f.write(TRP_ARCH_T_VOLCANIUM_NUC)
    #######################################################################

    # TODO Add proper aa data for val_euk_M_domestica into source data set
    directory = os.path.join(OUTPUT_DIR, 'Eukaryotes', 'M_domestica')
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename = os.path.join(directory, 'val_euk_M_domestica.aa')
    with open(filename, 'w') as f:
        f.write(VAL_EUK_M_DOMESTICA_AA)
    #######################################################################

    # TODO Add proper nuc data for asn_euk_G_lamblia into source data set
    directory = os.path.join(OUTPUT_DIR, 'Eukaryotes', 'G_lamblia')
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename = os.path.join(directory, 'asn_euk_G_lamblia.nuc')
    with open(filename, 'w') as f:
        f.write(GIARDIA_LAMBLIA_ASN_NUC)
    #######################################################################


LEU_ARCH_T_VOLCANIUM_AA = '''>gi|48428561
MDYSRCCGSTNSRVIYMNIDEKWQNAWERDHVFEPKIDERKKFMITVPWPYTNGSLHVGHGRTYTLGDII
ARYKRSRNYNVLFPMGFHQSGTPILAFSERIRAGDASTIALYRSYLSEYGEKDIDGWLEKFKDPRNIADY
FSNAIINDFKHLGYSIDWTRKFTSADEFYQNVVKWQFHKLNEKGLIKQDKYPILYSIDDDNAVGEDDIKD
GDTDKVSVEEYTAVFFESNSYSLIAASLRPETLFGVTNIWINPTGEYVKIKIGDKIAVVSKEAVDKLKYQ
RNDVSVIGPISAESIQRKKFTTPFGKEVPVYKADFVDTDNGTGVVYSVPSHSVYDFVYYRRKKSGQTPVV
IEAPLKMPEVEIKFDLNSKEGLDEATKELYKSEFYYGKLVNSGEYTGLTVRDAREKIKKDLIGSGKAIIF
YETSRKAVTRGGSKVIVAVLPDQWFIDYSADWLKKLSHDMLNRMMIYPEMYRNVMNDAIDWLKERPCARR
RGLGTKLPFDDRWVIESLSDSTIYPAVYTTSIQMRKLYENGKLDENAIERIFDGGEVQNDEERTARNEFS
YWYPVDIRLTAVPHISNHLSFYVMNHAAIFPPEKWPSGLIISGLVVSNGAKISKSKGNVVSLLEITKKYS
ADIYRLYVAVQADVSSTMDWNENDLSNIVRRFNEFKTIMDSFKPDTSELNFEETWFVSRFAERLKQFMDQ
MDGFQIRDAYINIFYGTLNDLKYAVNRGASQNRSLASIIADWLRALMPVISHHAEEYWHRYVSNTYVSIN
PFDDNFAEKYERLAKVYGLSTSEFYQVMDYVEHIIQDINNIISVTGIEPKSVEITVANEDVIKASREFLS
NSVSERSKRYLQYLAKRRKDIVVYPFNEIDILRRNSSYISRQVKADVSINTGDIINGKIAVPGKPVIHIT
'''

TRP_ARCH_T_VOLCANIUM_AA = '''>sp|Q978Y8|SYW_THEVO Tryptophan--tRNA ligase OS=Thermoplasma volcanium (strain ATCC 51530 / DSM 4299 / JCM 9571 / NBRC 15438 / GSS1) OX=273116 GN=trpS PE=3 SV=1
MINPWSSSDFFDYERLKKEFGISDQSDNIDHFLFRRKVILGQRGFEYIKYAIDNKIKFNV
MTGLMPSGEMHLGNKSAIDQVIYFQKLGGSVSIAVADLESYSTRGIPLDKAREIAIEKYI
LNYIAMGLQPCEIYFQSKNKDVQFLSYILGNWTNMNELKALYGFTDSNDILHINAPLIQA
ADVLHTQLNNYGGPAPTVVPVGFDQDPHIRLMRDLAKRMRIFNVFYDGGITVSIKGKGDS
TMPVDQAYEYLSKRFSEVTKDYEYRVVKAKDGKEEDIVRTDIDLAKIGSEFNVFSFIPPS
ATYQKLMKGLKGGKMSSSVPDSLISMNDDVEEAKRKIMRALTGGRDTEEEQRKLGGEPEK
CPVFDLYNYEIDDDKYVNEVFEECKSGKRMCGYCKREIADKMSIFLKDIKEKREIAREKL
SLYIHE
'''

TRP_ARCH_T_VOLCANIUM_NUC = '''
ATGATAAATCCGTGGTCTTCGAGCGACTTCTTTGATTATGAGAGGTTGAAGAAGGAATTTGGAATTTCAG
ATCAATCCGATAATATTGATCATTTCTTGTTTCGACGGAAGGTGATTCTCGGTCAGCGTGGCTTTGAATA
CATTAAGTATGCAATTGATAATAAGATCAAGTTCAACGTAATGACCGGGTTAATGCCGTCAGGAGAAATG
CATCTGGGAAATAAGAGTGCTATCGATCAAGTTATTTACTTTCAGAAGCTTGGAGGAAGCGTATCTATAG
CGGTAGCTGATCTTGAATCGTACTCTACTAGGGGTATACCGCTCGATAAAGCAAGGGAAATAGCTATTGA
GAAATACATACTTAACTACATCGCTATGGGCCTTCAACCTTGTGAGATCTATTTTCAATCTAAAAACAAA
GATGTTCAGTTTCTATCCTATATCCTTGGAAATTGGACTAATATGAATGAGCTCAAAGCTCTGTACGGTT
TTACTGATTCAAATGATATCTTACACATCAATGCGCCTCTCATCCAGGCAGCTGACGTTCTCCATACACA
GCTGAACAATTATGGCGGGCCAGCACCCACTGTAGTTCCTGTAGGTTTTGATCAGGATCCACATATAAGG
CTCATGAGGGATCTAGCTAAAAGAATGAGGATTTTCAATGTCTTCTATGATGGTGGGATTACAGTTTCAA
TAAAAGGGAAGGGAGATTCGACGATGCCTGTAGATCAGGCTTATGAATACCTTTCAAAGAGATTTTCTGA
AGTGACAAAAGACTATGAATACAGAGTGGTAAAGGCCAAGGATGGCAAGGAAGAAGATATTGTTAGAACA
GACATCGATCTTGCTAAGATCGGATCAGAGTTTAATGTATTTTCATTCATTCCTCCATCTGCAACCTACC
AGAAACTTATGAAGGGTTTAAAAGGAGGAAAGATGTCTTCATCGGTTCCAGATTCTCTTATTTCGATGAA
TGATGATGTGGAGGAAGCTAAAAGGAAAATAATGCGTGCGCTAACTGGAGGCAGGGATACTGAAGAAGAA
CAGAGAAAACTGGGAGGTGAACCTGAAAAATGCCCAGTCTTCGATCTATACAACTATGAAATAGACGATG
ACAAATACGTAAACGAAGTATTTGAAGAGTGCAAATCGGGAAAGAGGATGTGCGGATACTGCAAGAGGGA
AATAGCCGATAAAATGTCCATATTTTTAAAAGATATAAAGGAAAAAAGGGAAATTGCAAGAGAAAAATTA
TCACTTTATATTCATGAATAG
'''

VAL_EUK_M_DOMESTICA_AA = '''MSILYVSPHPDAFPSLRALIAARYGESGPGPGWGGPPPRVCLQP
PPTSGSCIPPPRLPVLEQGPGGLRVWGAAAVAQLLWPAGMGGPGGSRGATLVQQWVSY
ADGELVPAACGATLPALGLRNPTQDPQAALGALGRALGPLEERLRLHTYLAGEAPTLA
DLAAVTALLLPFRYVLDPPARWAWGNVSRWFMTCVQQPEFRAVLGEVSLCSGIRPIPQ
QPGTEASGPPKTAAQLKKEAKKREKLEKFQQKQKNQLHQPPPGEKKLKIEKKEKRDPG
VITYDIPTPPGEKKDVSGPMPDSYSPQYVEAAWYSWWENKGFFKPEYGRASLTEPNPR
GTFMMCIPPPNVTGSLHLGHALTNAIQDSLTRWHRMRGETTLWNPGCDHAGIATQVVV
EKKLWRERGMSRHQLGREAFLREVWKWKNEKGDRIYHQLKKLGGSLDWDRACFTMDPK
LSAAVTEAFVRLHNDGVIYRSTRLVNWSCSLNSAISDIEVDKKELSGRTLLSVPGYEE
KVEFGVIVSFAYKIEDSESNEEVVVATTRIETMLGDVAVAVHPNDPRYQHLRGKSVMH
PFLLRSLPIIFDEFVDMEFGTGAVKITPAHDQNDYEVGQRHKLEAVSIMDHRGNLINV
PPPFLGLPRFEARKAVLAALKDKGLFREVKDNPMVVPLCNRSKDVVEPLLKPQWYVRC
GEMAQAASAAVTRGDLKILPEVHQKIWHIWMDNIRDWCISRQLWWGHRIPAYFVTVND
PAVPPGEDPDGRYWVSGRNEEEAREKAAKEFGVPPDKISLSQDEDVLDTWFSSGLFPF
SILGWPNQTEDLSIFYPGTLLETGHDILFFWVARMVMLGLKLTGKLPFKEVYLHALVR
DAHGRKMSKSLGNVIDPLDVISGLSLQGLHDQLLNSNLDPSEMEKAKEGQKADFPNGI
PECGTDALRFGLCAYTSQGRDINLDVNRILGYRHFCNKLWNATKFALRALGDGFVPSP
TPQASSQESLADRWIRSRLSEAVGLSHQGFQAYDFPTITTAQYSFWLYDLCNVYLECL
KPVLSGKDQVAAESARQTLYTCLDVGLRLLSPFMPFVTEELYQRLPRRGPQAPPSLCV
TPYPEPDELSWKDPEAEAAFELALSITRAVRSLRADYNLTRSQPECFLEVADEATGTQ
ASAVSGYVQALSNTGAVNVLLPGSPAPQGCAVGLASDRCSVHLQLQGLVDPTRELAKL
RAKRGEAERQAQRLRERRAVPDYATKVPSQVQESEEAKLQQTEAELKKVDEAIALFEK
ML
'''

GIARDIA_LAMBLIA_ASN_NUC = '''>NW_002477099.1:331749-333395 Giardia lamblia ATCC 50803 SC_592, whole genome shotgun sequence
ATGGCAGACACTAAGAAAGCTCGTCAAGAGGACGACTGCATGGACGAGGATGTCAGCCATCCATACTACG
CCGCATACAAGAACGAAATACAGGCTGAGCTCGATGAGATAGCAACTATGCGGCAGGAGGATGGGAGTGA
GTTGCCAAATAAGCGTGTGAAGAAGCTCATCAAGACCATACATGCGAAATACAAGAAGCTCTACACAGCC
ACAAAGCAAGAAGAAAGGTCTGCCGATCAGCAGCGTCAGACAGAAGCACACCTTGCTGAGCTTCCCGAGC
GCATGGAGCCGTACACGGGTCTGGCAAAGCGTTTAAGAGTTAGGGACTTTCCTTATTACAAGGCAGTTGA
CAAGCCTGTTGTTGTGTGTGGTTGGTGCCACCGTGTGCGCGTCTCTTCTCCCAAACTTGCTTTCGTGGTA
CTTAGAGACGGGACCGGATACTGTCAACTGGTTCTTAACGAGGCCTGCATGGGTACTCGTAGGACTCAGA
CACTTCTGCGCACAGAGGCCGCTGTCAAAGCTGTAGGAATTCTTGTCGCTGATACGCGCGCACAAGGCGG
CTACGAAATTCAGTGTCTCTACTTTGACATTATTGGTCCGTCATCTGGAGAGTTTGAAACGCGGATCACC
CCAGAGTCTGGTCCCGATGCTCGTGCCAGAGAGCGCCATCTCATTCACCGAGGCGAACATGGGTCTGCAA
TTTTGAAAGCAAGGGGCCGGATTCTAAGGGCCTTCCGCAACCACTTCTACTTCAAAGGCTGGACAGAGGT
TACTCCACCAACAATAGTCAACACCGAATGCGAAGGAGGCTCTAGTCTGTTCAAAGTGGACTTTTATGGA
GAAGACGCATACCTTACCCAGTCCTCACAGCTTTATCTGGAGAGCGTCATTTCTGCCCTTGGGGATGTCT
ATTGCATCTTGCCTTCATTCAGAGCAGAGAAATCGAACACAAGGAGACACCTTTCCGAATTTACACATCT
TGAAACAGAGCACCCGTTCATTACGTTTTCGGAGCTCTTGGGGATAATTGAGGATTTCGTGATCTCTGTT
GTCACTGAGCTGAAAAATGACAAGCAATTCTATCCCCTTGTCCAGTTCCTTAATGGTGACAACTTTGACA
AGATTGCCTCAAGCTTGAAGAAACCCTTTAAACGTATGAGGCACAGTGAAGGAATAGACTGGCTGAATGA
TCACGGAATTACCAAGGAGGACGGTACAGCATTCACTTATGATGACGACATTCCAGAAGCCCAGGAGAGG
AAAATGATCGATGCCATAGGAGAACCGGTTTTCCTTACACACTTTCCAGCACATCTGAAGTCCTTTTACA
TGGCACGCTGCAGTGATGATCCCTCTTCCCCAGACTACCTTCTCACGGAAAGTGTCGACTTACTCGTTCC
CACTGTTGGCGAAATTGTCGGGGGGAGTATGCGTAAGTGGACGTACGAAGACACATACAAGGCCTGTGTT
GATGCAGGGCTCCCCATTAAGCGCTACTACTGGTATCTTGATTCACGCAAATTCGGTGGCTGTCCTCACG
GCGGAATGGGCCTGGGTCTGGAAAGGTTCATTTGCTGGTTACTAGGAATCTATACTGTTAGGGACACAGT
TCTTTACCCTCGTGCGCCGGGTCTCATTGAGCCGTGA'''


if __name__ == "__main__":
    load_missing_data()
    main()
