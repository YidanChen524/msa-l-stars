"""
helpers functions, classes and configurations
"""
import numpy as np
from collections import deque

# define gap and score matrix
gap = 5
score = np.array([[0, 5, 2, 5],
                  [5, 0, 5, 2],
                  [2, 5, 0, 5],
                  [5, 2, 5, 0]])

# define how nucleotides are mapped to indices
mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3,
           'a': 0, 'c': 1, 'g': 2, 't': 3,
           'N': 0, 'R': 0, 'S': 1}


def parse_fasta(filename):
    """helper functions for parsing fasta files"""
    with open(filename, 'r') as f:
        lines = f.readlines()
    names, seqs = [], []
    for i in range(len(lines)):
        if lines[i][0] == '>':
            names.append(lines[i][1:].strip('\n'))
            seq = ""
            while i < len(lines) - 1 and lines[i + 1][0] != '>':
                seq += lines[i + 1].strip('\n')
                i += 1
            seqs.append(seq)
    return names, seqs


def dynamic_table_2D(seq1, seq2, weight=1):
    """calculate the dynamic table between 2 sequences"""
    m, n = len(seq1) + 1, len(seq2) + 1
    t = np.zeros([m, n])
    for i in range(1, m):
        t[i, 0] = t[i - 1, 0] + gap
    for j in range(1, n):
        t[0, j] = t[0, j - 1] + gap
    for i in range(1, m):
        for j in range(1, n):
            v1 = t[i - 1, j] + gap
            v2 = t[i, j - 1] + gap
            v3 = t[i - 1, j - 1] + score[mapping[seq1[i - 1]], mapping[seq2[j - 1]]]
            t[i, j] = min(v1, v2, v3)
    return t * weight


def dynamic_table_3D(seq1, seq2, seq3, weight=1):
    """return the dynamic table of 3 sequences"""
    n1, n2, n3 = len(seq1) + 1, len(seq2) + 1, len(seq3) + 1
    t = np.zeros([n1, n2, n3])
    wg = weight * gap
    for i in range(1, n1):
        t[i, 0, 0] = t[i - 1, 0, 0] + 2 * wg
    for j in range(1, n2):
        t[0, j, 0] = t[0, j - 1, 0] + gap + wg
    for k in range(1, n3):
        t[0, 0, k] = t[0, 0, k - 1] + wg + gap
    for i in range(1, n1):
        for j in range(1, n2):
            v1 = t[i - 1, j, 0] + 2 * wg
            v2 = t[i, j - 1, 0] + wg + gap
            v3 = t[i - 1, j - 1, 0] + score[mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight + wg + gap
            t[i, j, 0] = min(v1, v2, v3)
    for i in range(1, n1):
        for k in range(1, n3):
            v1 = t[i - 1, 0, k] + 2 * wg
            v2 = t[i, 0, k - 1] + wg + gap
            v3 = t[i - 1, 0, k - 1] + score[mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight + wg + gap
            t[i, 0, k] = min(v1, v2, v3)
    for j in range(1, n2):
        for k in range(1, n3):
            v1 = t[0, j - 1, k] + wg + gap
            v2 = t[0, j, k - 1] + wg + gap
            v3 = t[0, j - 1, k - 1] + score[mapping[seq2[j - 1]], mapping[seq3[k - 1]]] + 2 * wg
            t[0, j, k] = min(v1, v2, v3)
    for i in range(1, n1):
        for j in range(1, n2):
            for k in range(1, n3):
                v1 = t[i - 1, j - 1, k - 1] + score[mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight \
                     + score[mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight \
                     + score[mapping[seq2[j - 1]], mapping[seq3[k - 1]]]
                v2 = t[i, j - 1, k - 1] + 2 * wg + score[mapping[seq2[j - 1]], mapping[seq3[k - 1]]]
                v3 = t[i - 1, j, k - 1] + wg + gap + score[
                    mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight
                v4 = t[i - 1, j - 1, k] + score[
                    mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight + wg + gap
                v5 = t[i, j, k - 1] + wg + gap
                v6 = t[i, j - 1, k] + wg + gap
                v7 = t[i - 1, j, k] + 2 * wg
                t[i, j, k] = min(v1, v2, v3, v4, v5, v6, v7)
    return t


def dynamic_table_5D_2l_star(seq0, seq1, seq2, seq3, seq4, weight=1):
    """
    return the dynamic table for 5 sequences, based on the graph configuration of (2l-1)-star.
    seq0 is the center string, (seq1, seq2) and (seq3, seq4) were in the same l-clique.
    """
    n0, n1, n2, n3, n4 = len(seq0) + 1, len(seq1) + 1, len(seq2) + 1, len(seq3) + 1, len(seq4) + 1
    t = np.zeros([n0, n1, n2, n3, n4])
    wg = weight * gap
    wg2, wg3, wg4 = wg * 2, wg * 3, wg * 4
    g2, g4 = gap * 2, gap * 4
    for i in range(1, n0):
        t[i, 0, 0, 0, 0] = t[i - 1, 0, 0, 0, 0] + wg4
    for j in range(1, n1):
        t[0, j, 0, 0, 0] = t[0, j - 1, 0, 0, 0] + g2 + wg
    for k in range(1, n2):
        t[0, 0, k, 0, 0] = t[0, 0, k - 1, 0, 0] + g2 + wg
    for l in range(1, n3):
        t[0, 0, 0, l, 0] = t[0, 0, 0, l - 1, 0] + g2 + wg
    for m in range(1, n4):
        t[0, 0, 0, 0, m] = t[0, 0, 0, 0, m - 1] + g2 + wg
    for i in range(1, n0):
        for j in range(1, n1):
            v1 = t[i - 1, j, 0, 0, 0] + wg4
            v2 = t[i, j - 1, 0, 0, 0] + wg + g2
            v3 = t[i - 1, j - 1, 0, 0, 0] + score[mapping[seq0[i - 1]], mapping[seq1[j - 1]]] * weight + wg2 + g2
            t[i, j, 0, 0, 0] = min(v1, v2, v3)
    for i in range(1, n0):
        for k in range(1, n2):
            v1 = t[i - 1, 0, k, 0, 0] + wg4
            v2 = t[i, 0, k - 1, 0, 0] + wg + g2
            v3 = t[i - 1, 0, k - 1, 0, 0] + score[mapping[seq0[i - 1]], mapping[seq2[k - 1]]] * weight + wg3 + g2
            t[i, 0, k, 0, 0] = min(v1, v2, v3)
    for i in range(1, n0):
        for l in range(1, n3):
            v1 = t[i - 1, 0, 0, l, 0] + wg4
            v2 = t[i, 0, 0, l - 1, 0] + wg + g2
            v3 = t[i - 1, 0, 0, l - 1, 0] + score[mapping[seq0[i - 1]], mapping[seq3[l - 1]]] * weight + wg3 + g2
            t[i, 0, 0, l, 0] = min(v1, v2, v3)
    for i in range(1, n0):
        for m in range(1, n4):
            v1 = t[i - 1, 0, 0, 0, m] + wg4
            v2 = t[i, 0, 0, 0, m - 1] + wg + g2
            v3 = t[i - 1, 0, 0, 0, m - 1] + score[mapping[seq0[i - 1]], mapping[seq4[m - 1]]] * weight + wg3 + g2
            t[i, 0, 0, 0, m] = min(v1, v2, v3)
    for j in range(1, n1):
        for k in range(1, n2):
            v1 = t[0, j - 1, k, 0, 0] + wg + g2
            v2 = t[0, j, k - 1, 0, 0] + wg + g2
            v3 = t[0, j - 1, k - 1, 0, 0] + wg2 + g4
            t[0, j, k, 0, 0] = min(v1, v2, v3)
    for j in range(1, n1):
        for l in range(1, n3):
            v1 = t[0, j - 1, 0, l, 0] + wg + g2
            v2 = t[0, j, 0, l - 1, 0] + wg + g2
            v3 = t[0, j - 1, 0, l - 1, 0] + score[mapping[seq1[j - 1]], mapping[seq3[l - 1]]] + wg2 + g2
            t[0, j, 0, l, 0] = min(v1, v2, v3)
    for j in range(1, n1):
        for m in range(1, n4):
            v1 = t[0, j - 1, 0, 0, m] + wg + g2
            v2 = t[0, j, 0, 0, m - 1] + wg + g2
            v3 = t[0, j - 1, 0, 0, m - 1] + score[mapping[seq1[j - 1]], mapping[seq4[m - 1]]] + wg2 + g2
            t[0, j, 0, 0, m] = min(v1, v2, v3)
    for k in range(1, n2):
        for l in range(1, n3):
            v1 = t[0, 0, k - 1, l, 0] + wg + g2
            v2 = t[0, 0, k, l - 1, 0] + wg + g2
            v3 = t[0, 0, k - 1, l - 1, 0] + score[mapping[seq2[k - 1]], mapping[seq3[l - 1]]] + wg2 + g2
            t[0, 0, k, l, 0] = min(v1, v2, v3)
    for k in range(1, n2):
        for m in range(1, n4):
            v1 = t[0, 0, k - 1, 0, m] + wg + g2
            v2 = t[0, 0, k, 0, m - 1] + wg + g2
            v3 = t[0, 0, k - 1, 0, m - 1] + score[mapping[seq2[k - 1]], mapping[seq4[m - 1]]] + wg2 + g2
            t[0, 0, k, 0, m] = min(v1, v2, v3)
    for l in range(1, n3):
        for m in range(1, n4):
            v1 = t[0, 0, 0, l - 1, m] + wg + g2
            v2 = t[0, 0, 0, l, m - 1] + wg + g2
            v3 = t[0, 0, 0, l - 1, m - 1] + wg2 + g4
            t[0, 0, 0, l, m] = min(v1, v2, v3)
    for i in range(1, n0):
        for j in range(1, n1):
            for k in range(1, n2):
                sij = score[mapping[seq0[i - 1]], mapping[seq1[j - 1]]]
                sik = score[mapping[seq0[i - 1]], mapping[seq2[k - 1]]]
                sjk = score[mapping[seq1[j - 1]], mapping[seq2[k - 1]]]
                v1 = t[i - 1, j - 1, k - 1, 0, 0] + (sij + sik) * weight + sjk + wg2 + g4
                v2 = t[i, j - 1, k - 1, 0, 0] + wg2 + sjk + g4
                v3 = t[i - 1, j, k - 1, 0, 0] + wg3 + g2 + sik * weight
                v4 = t[i - 1, j - 1, k, 0, 0] + sij * weight + wg3 + g2
                v5 = t[i, j, k - 1, 0, 0] + wg + g2
                v6 = t[i, j - 1, k, 0, 0] + wg + g2
                v7 = t[i - 1, j, k, 0, 0] + wg4
                t[i, j, k, 0, 0] = min(v1, v2, v3, v4, v5, v6, v7)
    for i in range(1, n0):
        for j in range(1, n1):
            for l in range(1, n3):
                sij = score[mapping[seq0[i - 1]], mapping[seq1[j-1]]]
                sil = score[mapping[seq0[i - 1]], mapping[seq3[l - 1]]]
                sjl = score[mapping[seq1[j - 1]], mapping[seq3[l - 1]]]
                v1 = t[i - 1, j - 1, 0, l - 1, 0] + (sij + sil) * weight + sjl + wg2 + g2
                v2 = t[i, j - 1, 0, l - 1, 0] + wg2 + sjl + g2
                v3 = t[i - 1, j, 0, l - 1, 0] + wg3 + g2 + sil * weight
                v4 = t[i - 1, j - 1, 0, l, 0] + sij * weight + wg3 + g2
                v5 = t[i, j, 0, l - 1, 0] + wg + g2
                v6 = t[i, j - 1, 0, l, 0] + wg + g2
                v7 = t[i - 1, j, 0, l, 0] + wg4
                t[i, j, 0, l, 0] = min(v1, v2, v3, v4, v5, v6, v7)
    for i in range(1, n0):
        for j in range(1, n1):
            for m in range(1, n4):
                sij = score[mapping[seq0[i - 1]], mapping[seq1[j-1]]]
                sim = score[mapping[seq0[i - 1]], mapping[seq4[m-1]]]
                sjm = score[mapping[seq1[j - 1]], mapping[seq4[m-1]]]
                v1 = t[i - 1, j - 1, 0, 0, m - 1] + (sij + sim) * weight + sjm + wg2 + g2
                v2 = t[i, j - 1, 0, 0, m - 1] + wg2 + sjm + g2
                v3 = t[i - 1, j, 0, 0, m - 1] + wg3 + g2 + sim * weight
                v4 = t[i - 1, j - 1, 0, 0, m] + sij * weight + wg3 + g2
                v5 = t[i, j, 0, 0, m - 1] + wg + g2
                v6 = t[i, j - 1, 0, 0, m] + wg + g2
                v7 = t[i - 1, j, 0, 0, m] + wg4
                t[i, j, 0, 0, m] = min(v1, v2, v3, v4, v5, v6, v7)
    for i in range(1, n0):
        for k in range(1, n2):
            for l in range(1, n3):
                sik = score[mapping[seq0[i - 1]], mapping[seq2[k-1]]]
                sil = score[mapping[seq0[i - 1]], mapping[seq3[l - 1]]]
                skl = score[mapping[seq2[k - 1]], mapping[seq3[l - 1]]]
                v1 = t[i - 1, 0, k - 1, l - 1, 0] + (sik + sil) * weight + skl + wg2 + g2
                v2 = t[i, 0, k - 1, l - 1, 0] + wg2 + skl + g2
                v3 = t[i - 1, 0, k, l - 1, 0] + wg3 + g2 + sil * weight
                v4 = t[i - 1, 0, k - 1, l, 0] + sik * weight + wg3 + g2
                v5 = t[i, 0, k, l - 1, 0] + wg + g2
                v6 = t[i, 0, k - 1, l, 0] + wg + g2
                v7 = t[i - 1, 0, k, l, 0] + wg4
                t[i, 0, k, l, 0] = min(v1, v2, v3, v4, v5, v6, v7)
    for i in range(1, n0):
        for k in range(1, n2):
            for m in range(1, n4):
                sik = score[mapping[seq0[i - 1]], mapping[seq2[k-1]]]
                sim = score[mapping[seq0[i - 1]], mapping[seq4[m-1]]]
                skm = score[mapping[seq2[k - 1]], mapping[seq4[m-1]]]
                v1 = t[i - 1, 0, k - 1, 0, m - 1] + (sik + sim) * weight + skm + wg2 + g2
                v2 = t[i, 0, k - 1, 0, m - 1] + wg2 + skm + g2
                v3 = t[i - 1, 0, k, 0, m - 1] + wg3 + g2 + sim * weight
                v4 = t[i - 1, 0, k - 1, 0, m] + sik * weight + wg3 + g2
                v5 = t[i, 0, k, 0, m - 1] + wg + g2
                v6 = t[i, 0, k - 1, 0, m] + wg + g2
                v7 = t[i - 1, 0, k, 0, m] + wg4
                t[i, 0, k, 0, m] = min(v1, v2, v3, v4, v5, v6, v7)
    for i in range(1, n0):
        for l in range(1, n3):
            for m in range(1, n4):
                sil = score[mapping[seq0[i - 1]], mapping[seq3[l - 1]]]
                sim = score[mapping[seq0[i - 1]], mapping[seq4[m-1]]]
                slm = score[mapping[seq3[l - 1]], mapping[seq4[m-1]]]
                v1 = t[i - 1, 0, 0, l - 1, m - 1] + (sil + sim) * weight + slm + wg2 + g4
                v2 = t[i, 0, 0, l - 1, m - 1] + wg2 + slm + g4
                v3 = t[i - 1, 0, 0, l, m - 1] + wg3 + g2 + sim * weight
                v4 = t[i - 1, 0, 0, l - 1, m] + sil * weight + wg3 + g2
                v5 = t[i, 0, 0, l, m - 1] + wg + g2
                v6 = t[i, 0, 0, l - 1, m] + wg + g2
                v7 = t[i - 1, 0, 0, l, m] + wg4
                t[i, 0, 0, l, m] = min(v1, v2, v3, v4, v5, v6, v7)
    for j in range(1, n1):
        for k in range(1, n2):
            for l in range(1, n3):
                sjl = score[mapping[seq1[j - 1]], mapping[seq3[l-1]]]
                skl = score[mapping[seq2[k - 1]], mapping[seq3[l-1]]]
                v1 = t[0, j - 1, k - 1, l - 1, 0] + sjl + skl + wg3 + g2
                v2 = t[0, j, k - 1, l - 1, 0] + wg2 + skl + g2
                v3 = t[0, j - 1, k, l - 1, 0] + wg2 + g2 + sjl
                v4 = t[0, j - 1, k - 1, l, 0] + wg2 + g4
                v5 = t[0, j, k, l - 1, 0] + wg + g2
                v6 = t[0, j, k - 1, l, 0] + wg + g2
                v7 = t[0, j - 1, k, l, 0] + wg + g2
                t[0, j, k, l, 0] = min(v1, v2, v3, v4, v5, v6, v7)
    for j in range(1, n1):
        for k in range(1, n2):
            for m in range(1, n4):
                sjm = score[mapping[seq1[j - 1]], mapping[seq4[m-1]]]
                skm = score[mapping[seq2[k - 1]], mapping[seq4[m-1]]]
                v1 = t[0, j - 1, k - 1, 0, m - 1] + sjm + skm + wg3 + g2
                v2 = t[0, j, k - 1, 0, m - 1] + wg2 + skm + g2
                v3 = t[0, j - 1, k, 0, m - 1] + wg2 + g2 + sjm
                v4 = t[0, j - 1, k - 1, 0, m] + wg2 + g4
                v5 = t[0, j, k, 0, m - 1] + wg + g2
                v6 = t[0, j, k - 1, 0, m] + wg + g2
                v7 = t[0, j - 1, k, 0, m] + wg + g2
                t[0, j, k, 0, m] = min(v1, v2, v3, v4, v5, v6, v7)
    for j in range(1, n1):
        for l in range(1, n3):
            for m in range(1, n4):
                sjl = score[mapping[seq1[j - 1]], mapping[seq3[l - 1]]]
                sjm = score[mapping[seq1[j - 1]], mapping[seq4[m-1]]]
                v1 = t[0, j - 1, 0, l - 1, m - 1] + sjl + sjm + wg3 + g2
                v2 = t[0, j, 0, l - 1, m - 1] + wg2 + g4
                v3 = t[0, j - 1, 0, l, m - 1] + wg2 + g2 + sjm
                v4 = t[0, j - 1, 0, l - 1, m] + wg2 + sjl + g2
                v5 = t[0, j, 0, l, m - 1] + wg + g2
                v6 = t[0, j, 0, l - 1, m] + wg + g2
                v7 = t[0, j - 1, 0, l, m] + wg + g2
                t[0, j, 0, l, m] = min(v1, v2, v3, v4, v5, v6, v7)
    for k in range(1, n2):
        for l in range(1, n3):
            for m in range(1, n4):
                skl = score[mapping[seq2[k - 1]], mapping[seq3[l - 1]]]
                skm = score[mapping[seq2[k - 1]], mapping[seq4[m-1]]]
                v1 = t[0, 0, k - 1, l - 1, m - 1] + skl + skm + wg3 + g2
                v2 = t[0, 0, k, l - 1, m - 1] + wg2 + g4
                v3 = t[0, 0, k - 1, l, m - 1] + wg2 + g2 + skm
                v4 = t[0, 0, k - 1, l - 1, m] + wg2 + skl + g2
                v5 = t[0, 0, k, l, m - 1] + wg + g2
                v6 = t[0, 0, k, l - 1, m] + wg + g2
                v7 = t[0, 0, k - 1, l, m] + wg + g2
                t[0, 0, k, l, m] = min(v1, v2, v3, v4, v5, v6, v7)
    for i in range(1, n0):
        for j in range(1, n1):
            for k in range(1, n2):
                for l in range(1, n3):
                    sij = score[mapping[seq0[i - 1]], mapping[seq1[j-1]]]
                    sik = score[mapping[seq0[i - 1]], mapping[seq2[k-1]]]
                    sil = score[mapping[seq0[i - 1]], mapping[seq3[l-1]]]
                    sjl = score[mapping[seq1[j - 1]], mapping[seq3[l-1]]]
                    skl = score[mapping[seq2[k - 1]], mapping[seq3[l-1]]]
                    v1 = t[i, j, k, l - 1, 0] + wg + g2
                    v2 = t[i, j, k - 1, l, 0] + wg + g2
                    v3 = t[i, j, k - 1, l - 1, 0] + wg2 + g2 + skl
                    v4 = t[i, j - 1, k, l, 0] + wg + g2
                    v5 = t[i, j - 1, k, l - 1, 0] + wg2 + sjl + g2
                    v6 = t[i, j - 1, k - 1, l, 0] + wg2 + g4
                    v7 = t[i, j - 1, k - 1, l - 1, 0] + wg3 + sjl + g2 + skl
                    v8 = t[i - 1, j, k, l, 0] + wg4
                    v9 = t[i - 1, j, k, l - 1, 0] + wg3 + weight * sil + g2
                    v10 = t[i - 1, j, k - 1, l, 0] + wg3 + weight * sik + g2
                    v11 = t[i - 1, j, k - 1, l - 1, 0] + wg2 + weight * (sik + sil) + g2 + skl
                    v12 = t[i - 1, j - 1, k, l, 0] + weight * sij + wg3 + g2
                    v13 = t[i - 1, j - 1, k, l - 1, 0] + weight * (sij + sil) + wg2 + g2 + sjl
                    v14 = t[i - 1, j - 1, k - 1, l, 0] + weight * (sij + sik) + wg2 + g4
                    v15 = t[i - 1, j - 1, k - 1, l - 1, 0] + weight * (sij + sik + sil) + wg + sjl + skl + g2
                    t[i, j, k, l, 0] = min(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15)
    for i in range(1, n0):
        for j in range(1, n1):
            for k in range(1, n2):
                for m in range(1, n4):
                    sij = score[mapping[seq0[i - 1]], mapping[seq1[j-1]]]
                    sik = score[mapping[seq0[i - 1]], mapping[seq2[k-1]]]
                    sim = score[mapping[seq0[i - 1]], mapping[seq4[m-1]]]
                    sjm = score[mapping[seq1[j - 1]], mapping[seq4[m-1]]]
                    skm = score[mapping[seq2[k - 1]], mapping[seq4[m-1]]]
                    v1 = t[i, j, k, 0, m - 1] + wg + g2
                    v2 = t[i, j, k - 1, 0, m] + wg + g2
                    v3 = t[i, j, k - 1, 0, m - 1] + wg2 + g2 + skm
                    v4 = t[i, j - 1, k, 0, m] + wg + g2
                    v5 = t[i, j - 1, k, 0, m - 1] + wg2 + sjm + g2
                    v6 = t[i, j - 1, k - 1, 0, m] + wg2 + g4
                    v7 = t[i, j - 1, k - 1, 0, m - 1] + wg3 + sjm + g2 + skm
                    v8 = t[i - 1, j, k, 0, m] + wg4
                    v9 = t[i - 1, j, k, 0, m - 1] + wg3 + weight * sim + g2
                    v10 = t[i - 1, j, k - 1, 0, m] + wg3 + weight * sik + g2
                    v11 = t[i - 1, j, k - 1, 0, m - 1] + wg2 + weight * (sik + sim) + g2 + skm
                    v12 = t[i - 1, j - 1, k, 0, m] + weight * sij + wg3 + g2
                    v13 = t[i - 1, j - 1, k, 0, m - 1] + weight * (sij + sim) + wg2 + sjm + g2
                    v14 = t[i - 1, j - 1, k - 1, 0, m] + weight * (sij + sik) + wg2 + g4
                    v15 = t[i - 1, j - 1, k - 1, 0, m - 1] + weight * (sij + sik + sim) + wg + sjm + skm + g2
                    t[i, j, k, 0, m] = min(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15)
    for i in range(1, n0):
        for j in range(1, n1):
            for l in range(1, n3):
                for m in range(1, n4):
                    sij = score[mapping[seq0[i - 1]], mapping[seq1[j-1]]]
                    sil = score[mapping[seq0[i - 1]], mapping[seq3[l-1]]]
                    sim = score[mapping[seq0[i - 1]], mapping[seq4[m-1]]]
                    sjl = score[mapping[seq1[j - 1]], mapping[seq3[l-1]]]
                    sjm = score[mapping[seq1[j - 1]], mapping[seq4[m-1]]]
                    v1 = t[i, j, 0, l, m - 1] + wg + g2
                    v2 = t[i, j, 0, l - 1, m] + wg + g2
                    v3 = t[i, j, 0, l - 1, m - 1] + wg2 + g4
                    v4 = t[i, j - 1, 0, l, m] + wg + g2
                    v5 = t[i, j - 1, 0, l, m - 1] + wg2 + sjm + g2
                    v6 = t[i, j - 1, 0, l - 1, m] + wg2 + g2 + sjl
                    v7 = t[i, j - 1, 0, l - 1, m - 1] + wg3 + sjl + g2 + sjm
                    v8 = t[i - 1, j, 0, l, m] + wg4
                    v9 = t[i - 1, j, 0, l, m - 1] + wg3 + weight * sim + g2
                    v10 = t[i - 1, j, 0, l - 1, m] + wg3 + weight * sil + g2
                    v11 = t[i - 1, j, 0, l - 1, m - 1] + wg2 + weight * (sil + sim) + g4
                    v12 = t[i - 1, j - 1, 0, l, m] + weight * sij + wg3 + g2
                    v13 = t[i - 1, j - 1, 0, l, m - 1] + weight * (sij + sim) + wg2 + sjm + g2
                    v14 = t[i - 1, j - 1, 0, l - 1, m] + weight * (sij + sil) + wg2 + g2 + sjl
                    v15 = t[i - 1, j - 1, 0, l - 1, m - 1] + weight * (sij + sil + sim) + wg + sjl + sjm + g2
                    t[i, j, 0, l, m] = min(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15)
    for i in range(1, n0):
        for k in range(1, n2):
            for l in range(1, n3):
                for m in range(1, n4):
                    sik = score[mapping[seq0[i - 1]], mapping[seq2[k-1]]]
                    sil = score[mapping[seq0[i - 1]], mapping[seq3[l-1]]]
                    sim = score[mapping[seq0[i - 1]], mapping[seq4[m-1]]]
                    skl = score[mapping[seq2[k - 1]], mapping[seq3[l-1]]]
                    skm = score[mapping[seq2[k - 1]], mapping[seq4[m-1]]]
                    v1 = t[i, 0, k, l, m - 1] + wg + g2
                    v2 = t[i, 0, k, l - 1, m] + wg + g2
                    v3 = t[i, 0, k, l - 1, m - 1] + wg2 + g4
                    v4 = t[i, 0, k - 1, l, m] + wg + g2
                    v5 = t[i, 0, k - 1, l, m - 1] + wg2 + skm + g2
                    v6 = t[i, 0, k - 1, l - 1, m] + wg2 + g2 + skl
                    v7 = t[i, 0, k - 1, l - 1, m - 1] + wg3 + skl + g2 + skm
                    v8 = t[i - 1, 0, k, l, m] + wg4
                    v9 = t[i - 1, 0, k, l, m - 1] + wg3 + weight * sim + g2
                    v10 = t[i - 1, 0, k, l - 1, m] + wg3 + weight * sil + g2
                    v11 = t[i - 1, 0, k, l - 1, m - 1] + wg2 + weight * (sil + sim) + g4
                    v12 = t[i - 1, 0, k - 1, l, m] + weight * sik + wg3 + g2
                    v13 = t[i - 1, 0, k - 1, l, m - 1] + weight * (sik + sim) + wg2 + skm + g2
                    v14 = t[i - 1, 0, k - 1, l - 1, m] + weight * (sik + sil) + wg2 + g2 + skl
                    v15 = t[i - 1, 0, k - 1, l - 1, m - 1] + weight * (sik + sil + sim) + wg + skl + skm + g2
                    t[i, 0, k, l, m] = min(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15)
    for j in range(1, n1):
        for k in range(1, n2):
            for l in range(1, n3):
                for m in range(1, n4):
                    sjl = score[mapping[seq1[j - 1]], mapping[seq3[l-1]]]
                    sjm = score[mapping[seq1[j - 1]], mapping[seq4[m-1]]]
                    skl = score[mapping[seq2[k - 1]], mapping[seq3[l-1]]]
                    skm = score[mapping[seq2[k - 1]], mapping[seq4[m-1]]]
                    v1 = t[0, j, k, l, m - 1] + wg + g2
                    v2 = t[0, j, k, l - 1, m] + wg + g2
                    v3 = t[0, j, k, l - 1, m - 1] + wg2 + g4
                    v4 = t[0, j, k - 1, l, m] + wg + g2
                    v5 = t[0, j, k - 1, l, m - 1] + wg2 + skm + g2
                    v6 = t[0, j, k - 1, l - 1, m] + wg2 + g2 + skl
                    v7 = t[0, j, k - 1, l - 1, m - 1] + wg3 + skl + g2 + skm
                    v8 = t[0, j - 1, k, l, m] + wg + g2
                    v9 = t[0, j - 1, k, l, m - 1] + wg2 + sjm + g2
                    v10 = t[0, j - 1, k, l - 1, m] + wg2 + sjl + g2
                    v11 = t[0, j - 1, k, l - 1, m - 1] + wg3 + sjl + sjm + g2
                    v12 = t[0, j - 1, k - 1, l, m] + wg2 + g4
                    v13 = t[0, j - 1, k - 1, l, m - 1] + wg3 + sjm + skm + g2
                    v14 = t[0, j - 1, k - 1, l - 1, m] + sjl + wg3 + g2 + skl
                    v15 = t[0, j - 1, k - 1, l - 1, m - 1] + wg4 + sjl + sjm + skl + skm
                    t[0, j, k, l, m] = min(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15)
    for i in range(1, n0):
        for j in range(1, n1):
            for k in range(1, n2):
                for l in range(1, n3):
                    for m in range(1, n4):
                        sij = score[mapping[seq0[i-1]], mapping[seq1[j-1]]]
                        sik = score[mapping[seq0[i-1]], mapping[seq2[k-1]]]
                        sil = score[mapping[seq0[i-1]], mapping[seq3[l-1]]]
                        sim = score[mapping[seq0[i-1]], mapping[seq4[m-1]]]
                        sjl = score[mapping[seq1[j-1]], mapping[seq3[l-1]]]
                        sjm = score[mapping[seq1[j-1]], mapping[seq4[m-1]]]
                        skl = score[mapping[seq2[k-1]], mapping[seq3[l-1]]]
                        skm = score[mapping[seq2[k-1]], mapping[seq4[m-1]]]
                        v1 = t[i, j, k, l, m - 1] + wg + g2
                        v2 = t[i, j, k, l - 1, m] + wg + g2
                        v3 = t[i, j, k, l - 1, m - 1] + wg2 + g4
                        v4 = t[i, j, k - 1, l, m] + wg + g2
                        v5 = t[i, j, k - 1, l, m - 1] + wg2 + g2 + skm
                        v6 = t[i, j, k - 1, l - 1, m] + wg2 + g2 + skl
                        v7 = t[i, j, k - 1, l - 1, m - 1] + wg3 + g2 + skl + skm
                        v8 = t[i, j - 1, k, l, m] + wg + g2
                        v9 = t[i, j - 1, k, l, m - 1] + wg2 + g2 + sjm
                        v10 = t[i, j - 1, k, l - 1, m] + wg2 + sjl + g2
                        v11 = t[i, j - 1, k, l - 1, m - 1] + wg3
                        v12 = t[i, j - 1, k - 1, l, m] + wg2 + g4
                        v13 = t[i, j - 1, k - 1, l, m - 1] + wg3 + g2 + sjm + skm
                        v14 = t[i, j - 1, k - 1, l - 1, m] + wg3 + sjl + g2 + skl
                        v15 = t[i, j - 1, k - 1, l - 1, m - 1] + wg4 + sjl + sjm + skl + skm
                        v16 = t[i - 1, j, k, l, m] + wg4
                        v17 = t[i - 1, j, k, l, m - 1] + wg3 + sim + g2
                        v18 = t[i - 1, j, k, l - 1, m] + wg3 + weight * sil + g2
                        v19 = t[i - 1, j, k, l - 1, m - 1] + wg2 + weight * (sil + sim) + g4
                        v20 = t[i - 1, j, k - 1, l, m] + wg3 + weight * sik + g2
                        v21 = t[i - 1, j, k - 1, l, m - 1] + wg2 + weight * (sik + sim) + g2 + skm
                        v22 = t[i - 1, j, k - 1, l - 1, m] + wg2 + weight * (sik + sil) + g2 + skl
                        v23 = t[i - 1, j, k - 1, l - 1, m - 1] + wg + weight * (sik + sil + sim) + g2 + skl + skm
                        v24 = t[i - 1, j - 1, k, l, m] + weight * sij + wg3 + g2
                        v25 = t[i - 1, j - 1, k, l, m - 1] + weight * (sij + sim) + wg2 + g2 + sjm
                        v26 = t[i - 1, j - 1, k, l - 1, m] + weight * (sij + sil) + wg2 + sjl + g2
                        v27 = t[i - 1, j - 1, k, l - 1, m - 1] + weight * (sij + sil + sim) + wg + sjl + sjm + g2
                        v28 = t[i - 1, j - 1, k - 1, l, m] + weight * (sij + sik) + wg2 + g4
                        v29 = t[i - 1, j - 1, k - 1, l, m - 1] + weight * (sij + sik + sim) + wg + g2 + sjm + skm
                        v30 = t[i - 1, j - 1, k - 1, l - 1, m] + weight * (sij + sik + sil) + wg + sjl + g2 + skl
                        v31 = t[i - 1, j - 1, k - 1, l - 1, m - 1] + weight * (sij + sik + sil + sim) + sjl + sjm + skl + skm
                        t[i, j, k, l, m] = min(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16, v17, v18, v19, v20, v21, v22, v23, v24, v25, v26, v27, v28, v29, v30, v31)
    return t


def pairwise_alignment(seq1, seq2, weight=1):
    """return the optimal alignment between 2 sequences"""
    # fill out the dynamic table t
    t = dynamic_table_2D(seq1, seq2, weight)
    weighted_gap = weight * gap
    # compute an alignment
    i, j = len(seq1), len(seq2)
    a1, a2 = deque(), deque()
    while i > 0 or j > 0:
        v = t[i, j]
        if i > 0 and j > 0 and v == t[i - 1, j - 1] + weight * score[mapping[seq1[i - 1]], mapping[seq2[j - 1]]]:
            a1.appendleft(seq1[i - 1])
            a2.appendleft(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and v == t[i - 1, j] + weighted_gap:
            a1.appendleft(seq1[i - 1])
            a2.appendleft('-')
            i -= 1
        elif j > 0 and v == t[i, j - 1] + weighted_gap:
            a2.appendleft(seq2[j - 1])
            a1.appendleft('-')
            j -= 1
    return a1, a2


def three_exact_alignment(seq1, seq2, seq3, weight=1):
    """return the optimal alignment between 3 sequences"""
    # fill out the dynamic table t
    t = dynamic_table_3D(seq1, seq2, seq3, weight)
    weighted_gap = weight * gap
    # compute an alignment
    i, j, k = len(seq1), len(seq2), len(seq3)
    a1, a2, a3 = deque(), deque(), deque()
    while i > 0 or j > 0 or k > 0:
        v = t[i, j, k]
        if i > 0 and j > 0 and k > 0 and \
                v == t[i - 1, j - 1, k - 1] + score[mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight \
                + score[mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight \
                + score[mapping[seq2[j - 1]], mapping[seq3[k - 1]]]:
            a1.appendleft(seq1[i - 1])
            a2.appendleft(seq2[j - 1])
            a3.appendleft(seq3[k - 1])
            i, j, k = i - 1, j - 1, k - 1
        elif j > 0 and k > 0 and v == t[i, j - 1, k - 1] + 2 * weighted_gap + \
                score[mapping[seq2[j - 1]], mapping[seq3[k - 1]]]:
            a1.appendleft('-')
            a2.appendleft(seq2[j - 1])
            a3.appendleft(seq3[k - 1])
            j, k = j - 1, k - 1
        elif i > 0 and k > 0 and \
                v == t[i - 1, j, k - 1] + weighted_gap + gap + score[
            mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight:
            a1.appendleft(seq1[i - 1])
            a2.appendleft('-')
            a3.appendleft(seq3[k - 1])
            i, k = i - 1, k - 1
        elif i > 0 and j > 0 and \
                v == t[i - 1, j - 1, k] + score[
            mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight + weighted_gap + gap:
            a1.appendleft(seq1[i - 1])
            a2.appendleft(seq2[j - 1])
            a3.appendleft('-')
            i, j = i - 1, j - 1
        elif k > 0 and v == t[i, j, k - 1] + weighted_gap + gap:
            a1.appendleft('-')
            a2.appendleft('-')
            a3.appendleft(seq3[k - 1])
            k -= 1
        elif j > 0 and v == t[i, j - 1, k] + weighted_gap + gap:
            a1.appendleft('-')
            a2.appendleft(seq2[j - 1])
            a3.appendleft('-')
            j -= 1
        elif i > 0 and v == t[i - 1, j, k] + 2 * weighted_gap:
            a1.appendleft(seq1[i - 1])
            a2.appendleft('-')
            a3.appendleft('-')
            i -= 1
    return a1, a2, a3


def five_exact_alignment_2l_star(seq0, seq1, seq2, seq3, seq4, weight):
    """return the optimal alignment for 5 sequences, based on the graph configuration of (2l-1)-star"""
    # fill out the dynamic table
    t = dynamic_table_5D_2l_star(seq0, seq1, seq2, seq3, seq4, weight)
    wg = weight * gap
    wg2, wg3, wg4 = wg * 2, wg * 3, wg * 4
    g2, g4 = gap * 2, gap * 4
    # compute an alignment
    a0, a1, a2, a3, a4 = deque(), deque(), deque(), deque(), deque()
    i, j, k, l, m = len(seq0), len(seq1), len(seq2), len(seq3), len(seq4)
    while i > 0 or j > 0 or k > 0 or l > 0 or m > 0:
        v = t[i, j, k, l, m]
        sij = score[mapping[seq0[i-1]], mapping[seq1[j-1]]]
        sik = score[mapping[seq0[i-1]], mapping[seq2[k-1]]]
        sil = score[mapping[seq0[i-1]], mapping[seq3[l-1]]]
        sim = score[mapping[seq0[i-1]], mapping[seq4[m-1]]]
        sjl = score[mapping[seq1[j-1]], mapping[seq3[l-1]]]
        sjm = score[mapping[seq1[j-1]], mapping[seq4[m-1]]]
        skl = score[mapping[seq2[k-1]], mapping[seq3[l-1]]]
        skm = score[mapping[seq2[k-1]], mapping[seq4[m-1]]]
        if m > 0 and v == t[i, j, k, l, m - 1] + wg + g2:
            a0.appendleft('-')
            a1.appendleft('-')
            a2.appendleft('-')
            a3.appendleft('-')
            a4.appendleft(seq4[m-1])
            m -= 1
        elif l > 0 and v == t[i, j, k, l - 1, m] + wg + g2:
            a0.appendleft('-')
            a1.appendleft('-')
            a2.appendleft('-')
            a3.appendleft(seq3[l-1])
            a4.appendleft('-')
            l -= 1
        elif l > 0 and m > 0 and v == t[i, j, k, l - 1, m - 1] + wg2 + g4:
            a0.appendleft('-')
            a1.appendleft('-')
            a2.appendleft('-')
            a3.appendleft(seq3[l - 1])
            a4.appendleft(seq4[m - 1])
            l -= 1
            m -= 1
        elif k > 0 and v == t[i, j, k - 1, l, m] + wg + g2:
            a0.appendleft('-')
            a1.appendleft('-')
            a2.appendleft(seq2[k-1])
            a3.appendleft('-')
            a4.appendleft('-')
            k -= 1
        elif k > 0 and m > 0 and v == t[i, j, k - 1, l, m - 1] + wg2 + g2 + skm:
            a0.appendleft('-')
            a1.appendleft('-')
            a2.appendleft(seq2[k-1])
            a3.appendleft('-')
            a4.appendleft(seq4[m-1])
            k -= 1
            m -= 1
        elif k > 0 and l > 0 and v == t[i, j, k - 1, l - 1, m] + wg2 + g2 + skl:
            a0.appendleft('-')
            a1.appendleft('-')
            a2.appendleft(seq2[k - 1])
            a3.appendleft(seq3[l - 1])
            a4.appendleft('-')
            k -= 1
            l -= 1
        elif k > 0 and l > 0 and m > 0 and v == t[i, j, k - 1, l - 1, m - 1] + wg3 + g2 + skl + skm:
            a0.appendleft('-')
            a1.appendleft('-')
            a2.appendleft(seq2[k - 1])
            a3.appendleft(seq3[l - 1])
            a4.appendleft(seq4[m - 1])
            k -= 1
            l -= 1
            m -= 1
        elif j > 0 and v == t[i, j - 1, k, l, m] + wg + g2:
            a0.appendleft('-')
            a1.appendleft(seq1[j - 1])
            a2.appendleft('-')
            a3.appendleft('-')
            a4.appendleft('-')
            j -= 1
        elif j > 0 and m > 0 and v == t[i, j - 1, k, l, m - 1] + wg2 + g2 + sjm:
            a0.appendleft('-')
            a1.appendleft(seq1[j - 1])
            a2.appendleft('-')
            a3.appendleft('-')
            a4.appendleft(seq4[m - 1])
            j -= 1
            m -= 1
        elif j > 0 and l > 0 and v == t[i, j - 1, k, l - 1, m] + wg2 + sjl + g2:
            a0.appendleft('-')
            a1.appendleft(seq1[j - 1])
            a2.appendleft('-')
            a3.appendleft(seq3[l - 1])
            a4.appendleft('-')
            j -= 1
            l -= 1
        elif j > 0 and l > 0 and m > 0 and v == t[i, j - 1, k, l - 1, m - 1] + wg3:
            a0.appendleft('-')
            a1.appendleft(seq1[j - 1])
            a2.appendleft('-')
            a3.appendleft(seq3[l - 1])
            a4.appendleft(seq4[m - 1])
            j -= 1
            l -= 1
            m -= 1
        elif j > 0 and k > 0 and v == t[i, j - 1, k - 1, l, m] + wg2 + g4:
            a0.appendleft('-')
            a1.appendleft(seq1[j - 1])
            a2.appendleft(seq2[k - 1])
            a3.appendleft('-')
            a4.appendleft('-')
            j -= 1
            k -= 1
        elif j > 0 and k > 0 and m > 0 and v == t[i, j - 1, k - 1, l, m - 1] + wg3 + g2 + sjm + skm:
            a0.appendleft('-')
            a1.appendleft(seq1[j - 1])
            a2.appendleft(seq2[k - 1])
            a3.appendleft('-')
            a4.appendleft(seq4[m - 1])
            j -= 1
            k -= 1
            m -= 1
        elif j > 0 and k > 0 and l > 0 and v == t[i, j - 1, k - 1, l - 1, m] + wg3 + sjl + g2 + skl:
            a0.appendleft('-')
            a1.appendleft(seq1[j - 1])
            a2.appendleft(seq2[k - 1])
            a3.appendleft(seq3[l - 1])
            a4.appendleft('-')
            j -= 1
            k -= 1
            l -= 1
        elif j > 0 and k > 0 and l > 0 and m > 0 and v == t[i, j - 1, k - 1, l - 1, m - 1] + wg4 + sjl + sjm + skl + skm:
            a0.appendleft('-')
            a1.appendleft(seq1[j - 1])
            a2.appendleft(seq2[k - 1])
            a3.appendleft(seq3[l - 1])
            a4.appendleft(seq4[m - 1])
            j -= 1
            k -= 1
            l -= 1
            m -= 1
        elif i > 0 and v == t[i - 1, j, k, l, m] + wg4:
            a0.appendleft(seq0[i - 1])
            a1.appendleft('-')
            a2.appendleft('-')
            a3.appendleft('-')
            a4.appendleft('-')
            i -= 1
        elif i > 0 and m > 0 and v == t[i - 1, j, k, l, m - 1] + wg3 + sim + g2:
            a0.appendleft(seq0[i - 1])
            a1.appendleft('-')
            a2.appendleft('-')
            a3.appendleft('-')
            a4.appendleft(seq4[m - 1])
            i -= 1
            m -= 1
        elif i > 0 and l > 0 and v == t[i - 1, j, k, l - 1, m] + wg3 + weight * sil + g2:
            a0.appendleft(seq0[i - 1])
            a1.appendleft('-')
            a2.appendleft('-')
            a3.appendleft(seq3[l - 1])
            a4.appendleft('-')
            i -= 1
            l -= 1
        elif i > 0 and l > 0 and m > 0 and v == t[i - 1, j, k, l - 1, m - 1] + wg2 + weight * (sil + sim) + g4:
            a0.appendleft(seq0[i - 1])
            a1.appendleft('-')
            a2.appendleft('-')
            a3.appendleft(seq3[l - 1])
            a4.appendleft(seq4[m - 1])
            i -= 1
            l -= 1
            m -= 1
        elif i > 0 and k > 0 and v == t[i - 1, j, k - 1, l, m] + wg3 + weight * sik + g2:
            a0.appendleft(seq0[i - 1])
            a1.appendleft('-')
            a2.appendleft(seq2[k - 1])
            a3.appendleft('-')
            a4.appendleft('-')
            i -= 1
            k -= 1
        elif i > 0 and k > 0 and m > 0 and v == t[i - 1, j, k - 1, l, m - 1] + wg2 + weight * (sik + sim) + g2 + skm:
            a0.appendleft(seq0[i - 1])
            a1.appendleft('-')
            a2.appendleft(seq2[k - 1])
            a3.appendleft('-')
            a4.appendleft(seq4[m - 1])
            i -= 1
            k -= 1
            m -= 1
        elif i > 0 and k > 0 and l > 0 and v == t[i - 1, j, k - 1, l - 1, m] + wg2 + weight * (sik + sil) + g2 + skl:
            a0.appendleft(seq0[i - 1])
            a1.appendleft('-')
            a2.appendleft(seq2[k - 1])
            a3.appendleft(seq3[l - 1])
            a4.appendleft('-')
            i -= 1
            k -= 1
            l -= 1
        elif i > 0 and k > 0 and l > 0 and m > 0 and v == t[i - 1, j, k - 1, l - 1, m - 1] + wg + weight * (sik + sil + sim) + g2 + skl + skm:
            a0.appendleft(seq0[i - 1])
            a1.appendleft('-')
            a2.appendleft(seq2[k - 1])
            a3.appendleft(seq3[l - 1])
            a4.appendleft(seq4[m - 1])
            i -= 1
            k -= 1
            l -= 1
            m -= 1
        elif i > 0 and j > 0 and v == t[i - 1, j - 1, k, l, m] + weight * sij + wg3 + g2:
            a0.appendleft(seq0[i - 1])
            a1.appendleft(seq1[j - 1])
            a2.appendleft('-')
            a3.appendleft('-')
            a4.appendleft('-')
            i -= 1
            j -= 1
        elif i > 0 and j > 0 and m > 0 and v == t[i - 1, j - 1, k, l, m - 1] + weight * (sij + sim) + wg2 + g2 + sjm:
            a0.appendleft(seq0[i - 1])
            a1.appendleft(seq1[j - 1])
            a2.appendleft('-')
            a3.appendleft('-')
            a4.appendleft(seq4[m - 1])
            i -= 1
            j -= 1
            m -= 1
        elif i > 0 and j > 0 and l > 0 and v == t[i - 1, j - 1, k, l - 1, m] + weight * (sij + sil) + wg2 + sjl + g2:
            a0.appendleft(seq0[i - 1])
            a1.appendleft(seq1[j - 1])
            a2.appendleft('-')
            a3.appendleft(seq3[l - 1])
            a4.appendleft('-')
            i -= 1
            j -= 1
            l -= 1
        elif i > 0 and j > 0 and l > 0 and m > 0 and v == t[i - 1, j - 1, k, l - 1, m - 1] + weight * (sij + sil + sim) + wg + sjl + sjm + g2:
            a0.appendleft(seq0[i - 1])
            a1.appendleft(seq1[j - 1])
            a2.appendleft('-')
            a3.appendleft(seq3[l - 1])
            a4.appendleft(seq4[m - 1])
            i -= 1
            j -= 1
            l -= 1
            m -= 1
        elif i > 0 and j > 0 and k > 0 and v == t[i - 1, j - 1, k - 1, l, m] + weight * (sij + sik) + wg2 + g4:
            a0.appendleft(seq0[i - 1])
            a1.appendleft(seq1[j - 1])
            a2.appendleft(seq2[k - 1])
            a3.appendleft('-')
            a4.appendleft('-')
            i -= 1
            j -= 1
            k -= 1
        elif i > 0 and j > 0 and k > 0 and m > 0 and v == t[i - 1, j - 1, k - 1, l, m - 1] + weight * (sij + sik + sim) + wg + g2 + sjm + skm:
            a0.appendleft(seq0[i - 1])
            a1.appendleft(seq1[j - 1])
            a2.appendleft(seq2[k - 1])
            a3.appendleft('-')
            a4.appendleft(seq4[m - 1])
            i -= 1
            j -= 1
            k -= 1
            m -= 1
        elif i > 0 and j > 0 and k > 0 and l > 0 and v == t[i - 1, j - 1, k - 1, l - 1, m] + weight * (sij + sik + sil) + wg + sjl + g2 + skl:
            a0.appendleft(seq0[i - 1])
            a1.appendleft(seq1[j - 1])
            a2.appendleft(seq2[k - 1])
            a3.appendleft(seq3[l - 1])
            a4.appendleft('-')
            i -= 1
            j -= 1
            k -= 1
            l -= 1
        elif v == t[i - 1, j - 1, k - 1, l - 1, m - 1] + weight * (sij + sik + sil + sim) + sjl + sjm + skl + skm:
            a0.appendleft(seq0[i - 1])
            a1.appendleft(seq1[j - 1])
            a2.appendleft(seq2[k - 1])
            a3.appendleft(seq3[l - 1])
            a4.appendleft(seq4[m - 1])
            i -= 1
            j -= 1
            k -= 1
            l -= 1
            m -= 1
    return a0, a1, a2, a3, a4


def sp_score_clique(seqs, clique, k, l):
    """return the sp score of a clique"""
    if l == 2:
        return dynamic_table_2D(seqs[clique[0]], seqs[clique[1]], k - (l - 1))[-1, -1]
    if l == 3:
        return dynamic_table_3D(seqs[clique[0]], seqs[clique[1]], seqs[clique[2]], k - (l - 1))[-1, -1, -1]


def sp_score_clique_2l_star(seqs, clique, k, l):
    """return the sp scpre of a clique based on (2l-1)-star configuration"""
    if l == 2:
        return dynamic_table_3D(seqs[clique[0]], seqs[clique[1]], seqs[clique[2]], k - (l - 1) - 0.5)[-1, -1, -1]
    if l == 3:
        return dynamic_table_5D_2l_star(*[seqs[c] for c in clique], k - (l - 1) - 0.5)[-1, -1, -1, -1, -1]


def alignment_clique(seqs, clique, k, l):
    """return the optimal alignment of a clique"""
    if l == 2:
        return pairwise_alignment(seqs[clique[0]], seqs[clique[1]], k - (l - 1))
    if l == 3:
        return three_exact_alignment(seqs[clique[0]], seqs[clique[1]], seqs[clique[2]], k - (l - 1))


def alignment_clique_2l(seqs, clique, k, l):
    """return the optimal alignment of a 2l-1 clique"""
    if l == 2:
        return three_exact_alignment(seqs[clique[0]], seqs[clique[1]], seqs[clique[2]], k - (l - 1) - 0.5)
    if l == 3:
        return five_exact_alignment_2l_star(*[seqs[c] for c in clique], k - (l - 1) - 0.5)


def align_l_star(seqs, l_star, k, l):
    """given an l_star, return the optimal alignment of those sequences"""

    # a class that store a column of alignment
    class Column:
        def __init__(self):
            self.val = ['-'] * k
            self.next = None

    # initialize alignment as a linked list of Columns
    center = l_star[0][0]
    alignment = current = Column()
    for i in range(len(seqs[center])):
        current.next = Column()
        current = current.next
        current.val[center] = seqs[center][i]
    # merge cliques alignments
    for clique in l_star:
        a = alignment_clique(seqs, clique, k, l)
        current = alignment
        i = 0
        while i < len(a[0]):
            if not current.next:
                new_col = Column()
                for j in clique:
                    new_col.val[j] = a[clique.index(j)][i]
                new_col.next = current.next
                current.next = new_col
                current = current.next
                i += 1
            elif current.next.val[center] == a[0][i]:
                for j in clique:
                    current.next.val[j] = a[clique.index(j)][i]
                current = current.next
                i += 1
            elif current.next.val[center] == '-':
                current = current.next
            elif a[0][i] == '-':
                new_col = Column()
                for j in clique:
                    new_col.val[j] = a[clique.index(j)][i]
                new_col.next = current.next
                current.next = new_col
                current = current.next
                i += 1
    # return the alignment
    strings = [[] for _ in range(k)]
    current = alignment
    while current.next:
        for i in range(k):
            strings[i].append(current.next.val[i])
        current = current.next
    for i in range(k):
        strings[i] = "".join(strings[i])
    return strings


def align_2l_star(seqs, star, k, l):
    """given (2l-1)_star, return the optimal alignment of those sequences"""

    # a class that store a column of alignment
    class Column:
        def __init__(self):
            self.val = ['-'] * k
            self.next = None

    # initialize alignment as a linked list of Columns
    center = star[0][0]
    alignment = current = Column()
    for i in range(len(seqs[center])):
        current.next = Column()
        current = current.next
        current.val[center] = seqs[center][i]
    # merge cliques alignments
    for clique in star:
        a = alignment_clique_2l(seqs, clique, k, l)
        current = alignment
        i = 0
        while i < len(a[0]):
            if not current.next:
                new_col = Column()
                for j in clique:
                    new_col.val[j] = a[clique.index(j)][i]
                new_col.next = current.next
                current.next = new_col
                current = current.next
                i += 1
            elif current.next.val[center] == a[0][i]:
                for j in clique:
                    current.next.val[j] = a[clique.index(j)][i]
                current = current.next
                i += 1
            elif current.next.val[center] == '-':
                current = current.next
            elif a[0][i] == '-':
                new_col = Column()
                for j in clique:
                    new_col.val[j] = a[clique.index(j)][i]
                new_col.next = current.next
                current.next = new_col
                current = current.next
                i += 1
    # return the alignment
    strings = [[] for _ in range(k)]
    current = alignment
    while current.next:
        for i in range(k):
            strings[i].append(current.next.val[i])
        current = current.next
    for i in range(k):
        strings[i] = "".join(strings[i])
    return strings


if __name__ == "__main__":
    names, seqs = parse_fasta("test_seqs/testdata_7_seqs.txt")
    alignments = align_2l_star(seqs, [(2, 0, 3), (2, 6, 1), (2, 5, 4)], 7, 2)
    for alm in alignments:
        print(alm)
