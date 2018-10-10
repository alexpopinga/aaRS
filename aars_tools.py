"""Python tools for aaRS data

a.k.a. "Paul's tools"
"""
import os.path
import glob
import json
import xml.etree.ElementTree as ET
from urllib.parse import urlencode
from urllib.request import urlopen


try:
    from aars_algorithms_fast import levenshtein_distance_c as levenshtein_distance
    from aars_algorithms_fast import align_c as align
except ImportError:
    print('using slow algorithms')
    from aars_algorithms_slow import levenshtein_distance_py as levenshtein_distance
    from aars_algorithms_slow import align_py as align


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
    total = len(no_gaps)
    for i, no_gaps_key in enumerate(no_gaps.keys()):
        print('({:d}/{:d}) Aligning {}'.format(i+1, total, no_gaps_key))
        gaps_key, distance = find_match(
            no_gaps_key, no_gaps=no_gaps, gaps=gaps)
        domain, species, aa = parse_no_gaps_keys(no_gaps_key)
        if species == 'R_marinus':
            domain = 'arch'
        if species == 'R_rosetta':
            species = 'S_rosetta'
        path, path_cost = predict_path(domain, species, aa)
        aa_paths = [p for p in predict_amino_path(path, aa)
                    if "CUT" not in p and "TEST" not in p]

        # TODO What is this file?
        # Archaeans_Bacteria_Eukaryotes/Bacteria/Chroococcidiopsis thermalis/amino acid sequences/Cthermalis_asn_asp_aa
        # Need to fix and remove exception code.
        aa_paths = [p for p in aa_paths if 'Cthermalis_asn_asp_aa' not in p]
        ###############################################################################################################

        if len(aa_paths) > 1:
            raise RuntimeError('Too many paths found!\n{}'.format(aa_paths))
        try:
            aa_path = aa_paths[0]
        except IndexError:
            aa_path = None
        nuc_paths = predict_nucleotide_path(path, aa)
        if len(aa_paths) > 1:
            raise RuntimeError('Too many paths found!\n{}'.format(nuc_paths))
        try:
            nuc_path = nuc_paths[0]
        except IndexError:
            nuc_path = None
        aa_header, backup_id, aa_data = read_path_data(aa_path)
        if aa_data:
            aa_align = align(gaps[gaps_key], aa_data)
            aa_misalignments = sum([1 if c2 not in '-.*?' and c1 not in '*?' and c1 != c2 else 0
                                    for c1, c2 in zip(aa_data, aa_align)])
        else:
            aa_align = None
            aa_misalignments = None
        nuc_header, nuc_id, nuc_data = read_path_data(nuc_path)
        data[no_gaps_key] = {
            'gi': nuc_id if nuc_id else backup_id,
            'nuc_header': nuc_header,
            'aa-header': aa_header,
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
            'amino-acid-path': aa_path,
            'nucleotide-path': nuc_path,
            'amino-acid-data': aa_data,
            'amino-acid-aligned-gaps': aa_align,
            'amino-acid-misalignments': aa_misalignments,
            'nucleotide-data': nuc_data
        }
    with open('gap_data.json', 'w') as f:
        json.dump(data, f, indent=2)
    output_split_files()


def output_split_files():
    """Output files with summaries of aligned and unaligned"""
    with open('gap_data.json') as f:
        d = json.load(f)
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
        if not v['amino-acid-path'] and not v['nucleotide-path']:
            counts['missing data'] += 1
            missing_data[k] = {
                'name': k,
                'comment': "Missing amino acid and nucleotide sequence file"
            }
            continue
        if not v['nucleotide-path']:
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
        if misalignment == 0 and (estimated_nuc_len == actual_nuc_len or actual_nuc_len - estimated_nuc_len == 3):
            # From Peter Wills - 10 Oct 2018
            # "Many of these aaRS sequences are of the “catalytic domain” and some have
            # other domains attached to them, so you may have the first codon of the
            # following domain."
            #
            # Given this, if the nucleotide sequences has one extra codon, we ignore it.
            counts['perfectly aligned to our nucleotide data'] += 1
            nuc_align = ''.join(['---' if v['amino-acid-aligned-gaps'][i] in '-.'
                                 else v['nucleotide-data'][i*3:(i+1)*3]
                                 for i in range(len(v['amino-acid-aligned-gaps']))])
            perfect[k] = {
                'name': k,
                'regions': v['gaps-value'],
                'amino-acid-file': v['amino-acid-path'],
                'nucleotide-file': v['nucleotide-path'],
                'aa-num': ''.join(['{0:^3d}'.format(n) for n in range(len(v['amino-acid-data']))]),
                'aa-seq': ' ' + '  '.join(list(v['amino-acid-data'])) + ' ',
                'aa-ali': ' ' + '  '.join(list(v['amino-acid-aligned-gaps'])) + ' ',
                'nuc-al': nuc_align
            }
            continue
        # see if GenBank can do better
        check, gb_aa, gb_nuc = check_genbank(v['gi'])
        if check:
            aa_align = align(v['gaps-value'], gb_aa)
            aa_misalignments = sum([1 if c2 not in '-.*?' and c1 not in '*?' and c1 != c2 else 0
                                    for c1, c2 in zip(gb_aa, aa_align)])
            if aa_misalignments == 0:
                counts['perfectly aligned to GenBank nucleotide data'] += 1
                nuc_align = ''.join(['---' if v['amino-acid-aligned-gaps'][i] in '-.'
                                     else v['nucleotide-data'][i*3:(i+1)*3]
                                     for i in range(len(v['amino-acid-aligned-gaps']))])
                genbank[k] = {
                    'name': k,
                    'regions': v['gaps-value'],
                    'amino-acid-file': 'GenBank ' + v['gi'],
                    'nucleotide-file': 'GenBank ' + v['gi'],
                    'aa-num': ''.join(['{0:^3d}'.format(n) for n in range(len(gb_aa))]),
                    'aa-seq': ' ' + '  '.join(list(gb_aa)) + ' ',
                    'aa-ali': ' ' + '  '.join(list(aa_align)) + ' ',
                    'nuc-al': gb_nuc
                }
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
                'genbank-id': v['gi'],
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
    for k, v in counts.items():
        print('{}: {}'.format(k, v))
    with open('gap_data_perfect.txt', 'w') as p:
        json.dump(perfect, p, indent=2)
    with open('gap_data_genbank.txt', 'w') as gb:
        json.dump(genbank, gb, indent=2)
    with open('gap_data_good.txt', 'w') as g:
        json.dump(good, g, indent=2)
    with open('gap_data_bad.txt', 'w') as b:
        json.dump(bad, b, indent=2)
    with open('gap_data_missing.txt', 'w') as m:
        json.dump(missing_data, m, indent=2)


def predict_amino_path(path, aa):
    """Try to find the amino acid sequence file"""
    if path[-3:] == "ile" and aa == 'ile':  # M_mobile ends in 'ile' but may not be aaRS 'ile'
        return predict_amino_path(path, '_ile')
    paths = glob.glob(path + '/amino*/*{}_*'.format(aa))
    if not paths:
        paths = glob.glob(path + '/**/*{}_aa*'.format(aa), recursive=True)
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
    if path[-3:] == "ile" and aa == 'ile':  # M_modile ends in 'ile' but may not be aaRS 'ile'
        return predict_nucleotide_path(path, '_ile')
    paths = glob.glob(path + '/nuc*/*{}_*'.format(aa))
    if not paths:
        paths = glob.glob(path + '/**/*{}_nuc*'.format(aa), recursive=True)
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
    if path == None:
        return None, None, None
    dat = ''
    header = None
    gi = None
    with open(path) as path_p:
        for next_dat in path_p:
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
    nuc_data = get_genbank_nuc(gi)
    aa_data = get_genbank_aa(gi)
    nuc_len = len(nuc_data)
    aa_len = 3 * len(aa_data)
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


if __name__ == "__main__":
    make_no_gaps_to_gaps_file()
