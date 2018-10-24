"""Slow Python version of the algorithms"""
from itertools import product
import numpy as np


def levenshtein_distance_py(source, target):
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


def align_py(gapped_seq, full_seq):
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
