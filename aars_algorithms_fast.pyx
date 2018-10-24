import numpy as np
cimport numpy as np
cimport cython

def levenshtein_distance_c(source, target):
    source_bytes = source.encode('UTF-8')
    target_bytes = target.encode('UTF-8')

    cdef char* c_source = source_bytes
    cdef int c_source_length = len(source)
    cdef char* c_target = target_bytes
    cdef int c_target_length = len(target)
    cdef int[:, ::1] costs = np.empty((c_source_length+1, c_target_length+1), dtype='int32')

    return int(_calc_distance(c_source, c_source_length, c_target, c_target_length, costs))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef int _calc_distance(char* c_source, int c_source_length, char* c_target, int c_target_length, int[:, ::1] costs):
    cdef int i, j, s
    for i in range(c_source_length + 1):
        costs[i, 0] = i
    for j in range(c_target_length + 1):
        costs[0, j] = j
    for i in range(1, c_source_length + 1):
        for j in range(1, c_target_length + 1):
            if c_source[i-1] == c_target[j-1]:
                s = 0
            else:
                s = 1
            if costs[i-1, j] + 1 < costs[i, j-1] + 1 and costs[i-1, j] + 1 < costs[i-1, j-1] + s:
                costs[i, j] = costs[i-1, j] + 1
            elif costs[i, j-1] + 1 < costs[i-1, j-1] + s:
                costs[i, j] = costs[i, j-1] + 1
            else:
                costs[i, j] = costs[i-1, j-1] + s
    return costs[c_source_length, c_target_length]


def align_c(gapped_seq, full_seq):
    """align a gapped sequence to the full sequence"""
    regions = gapped_seq.strip('-') + '-'
    reg = regions.encode('UTF-8')
    seq = full_seq.encode('UTF-8')
    cdef char* reg_bytes = reg
    cdef int reg_len = len(regions)
    cdef char* seq_bytes = seq
    cdef int seq_len = len(full_seq)

    gaps = sum([1 if c == '-' else 0 for c in regions])
    cdef int matrix_side_1 = len(regions) - gaps + 1
    cdef int matrix_side_2 = (len(full_seq) - (len(regions) - gaps)) + 2

    if matrix_side_1 < 1 or matrix_side_2 < 1:
        return ''
    shape = (matrix_side_1, matrix_side_2)
    cdef int[:, ::1] costs = np.zeros(shape, dtype='int32')
    cdef char[:, ::1] path = np.zeros(shape, dtype='uint8')

    _align_c(reg_bytes, reg_len, seq_bytes, seq_len, costs, path, matrix_side_1, matrix_side_2)
    cdef int x, y
    x = matrix_side_1 - 1
    y = matrix_side_2 - 1
    alignment = []
    regions = [c for c in regions if c != '-']
    while x > 0 or y > 1:
        if path[x, y] == 1:  # take part
            alignment.append(regions.pop())
            x -= 1
        elif path[x, y] == 2:  # take '-' gap
            alignment.append('-')
            y -= 1
        else:  # take gap '.' gap
            alignment.append('.')
            y -= 1
    return ''.join(reversed(alignment))

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _align_c(char* reg, int reg_len, char* seq, int seq_len, int[:, ::1] costs, char[:, ::1] path, int s1, int s2):
    cdef int i, j, x, y
    cdef int gaps_cost, char_cost, take_gap, take_char
    cdef char c1, c2, c_star, c_question
    c_star = 42  # '*'
    c_question = 63
    for i in range(s1):
        costs[i, 0] = 999999999
        path[i, 0] = 0  # invalid
    for j in range(s2 - 1):
        costs[0, j + 1] = j
        path[0, j + 1] = 2  # take '-' gap
    i = 0
    for x in range(1, s1):
        c1 = reg[i]
        j = 0
        for y in range(1, s2):
            c2 = seq[j + x - 1]
            if reg[i + 1] == '-':
                gaps_cost = 1  # cost to extend a predefined gap
            elif path[x, y - 1] != 1:
                gaps_cost = 10  # cost to extend a new gap
            else:
                gaps_cost = 100  # cost to add a new gap
            take_gap = gaps_cost + costs[x, y - 1]
            if c1 == c2 or c1 == c_star or c1 == c_question or c2 == c_star or c_star == c_question:
                char_cost = 0  # cost to take a matching character
            else:
                char_cost = 10000  # cost to take a mismatched character
            take_char = char_cost + costs[x - 1, y]
            if take_char < take_gap:
                costs[x, y] = take_char
                path[x, y] = 1  # take character
            else:
                costs[x, y] = take_gap
                if gaps_cost == 1:
                    path[x, y] = 2  # take '-'
                else:
                    path[x, y] = 3  # take '.'
            j += 1
        i += 1
        while i < reg_len and reg[i] == '-':
            i += 1
