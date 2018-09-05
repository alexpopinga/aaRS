"""Read POA data and write FASTA file using all alignments"""
import os.path
import sys


def poa_to_fasta(filename):
    """Reads a POA file and writes a FASTA fail"""
    with open(filename) as data:
        pro_tag = False
        node_tag = False
        output_strings = []
        for line in data:
            if '<PRO>' in line:
                pro_tag = True
                continue
            elif pro_tag:
                if '</PRO>' in line:
                    pro_tag = False
                    continue
                output_strings.append(list('>' + line.split()[0]))
                output_strings.append([])
            elif '<LEN>' in line:
                length = int(line.split()[1])
                for i in range(1, len(output_strings), 2):
                    output_strings[i] = list('-' * length)
            elif '<NODE>' in line:
                node_tag = True
                continue
            elif node_tag:
                if '</NODE>' in line:
                    node_tag = False
                    break
                node = line.split()
                for i in range(1, len(node), 2):
                    index = int(node[0])
                    protein_num = int(node[i])
                    residue = node[i+1].split('.')[-1]
                    output_strings[2*protein_num+1][index] = residue
    with open(os.path.splitext(filename)[0] + '.fasta', 'w') as out:
        out.write('\n'.join([''.join(line) for line in output_strings]) + '\n')


if __name__ == "__main__":
    try:
        poa_to_fasta(sys.argv[1])
    except IndexError:
        print('Usage: python {} POA_FILENAME'.format(sys.argv[0]))
