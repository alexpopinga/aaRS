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