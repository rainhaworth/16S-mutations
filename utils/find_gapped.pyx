# relevant tutorial pages
# https://cython.readthedocs.io/en/latest/src/tutorial/numpy.html
# https://cython.readthedocs.io/en/latest/src/userguide/memoryviews.html
# https://cython.readthedocs.io/en/latest/src/tutorial/strings.html
import numpy as np
cimport numpy as cnp


cnp.import_array()

# input: sequence (string), kmer character (string) ndarray, offset (int) ndarray
def find_gapped(str seq, 
                cnp.ndarray kchars,
                cnp.ndarray[cnp.int64_t, ndim=2] koffs):
    cdef int seqlen = len(seq)
    cdef int kmer_count = koffs.shape[0]
    cdef int klen = koffs.shape[1]
    cdef int i, j, b, ksz, k_in_s
    cdef str subseq
    cdef bint found
    cdef cnp.ndarray counts_out = np.zeros([kmer_count], dtype=np.int64)

    for i in range(kmer_count):
        ksz = koffs[i][klen-1] + 1
        k_in_s = seqlen + 1 - ksz
        for j in range(k_in_s):
            subseq = seq[j:j+ksz]
            found = True
            for b in range(klen):
                if subseq[koffs[i][b]] != kchars[i][b]:
                    found = False
                    break
            if found:
                counts_out[i] = 1
                break

    return counts_out
