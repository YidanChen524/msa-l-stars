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
    weighted_gap = weight * gap
    for i in range(1, n1):
        t[i, 0, 0] = t[i - 1, 0, 0] + 2 * weighted_gap
    for j in range(1, n2):
        t[0, j, 0] = t[0, j - 1, 0] + gap + weighted_gap
    for k in range(1, n3):
        t[0, 0, k] = t[0, 0, k - 1] + weighted_gap + gap
    for i in range(1, n1):
        for j in range(1, n2):
            v1 = t[i - 1, j, 0] + 2 * weighted_gap
            v2 = t[i, j - 1, 0] + weighted_gap + gap
            v3 = t[i - 1, j - 1, 0] + score[mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight + weighted_gap + gap
            t[i, j, 0] = min(v1, v2, v3)
    for i in range(1, n1):
        for k in range(1, n3):
            v1 = t[i - 1, 0, k] + 2 * weighted_gap
            v2 = t[i, 0, k - 1] + weighted_gap + gap
            v3 = t[i - 1, 0, k - 1] + score[mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight + weighted_gap + gap
            t[i, 0, k] = min(v1, v2, v3)
    for j in range(1, n2):
        for k in range(1, n3):
            v1 = t[0, j - 1, k] + weighted_gap + gap
            v2 = t[0, j, k - 1] + weighted_gap + gap
            v3 = t[0, j - 1, k - 1] + score[mapping[seq2[j - 1]], mapping[seq3[k - 1]]] + 2 * weighted_gap
            t[0, j, k] = min(v1, v2, v3)
    for i in range(1, n1):
        for j in range(1, n2):
            for k in range(1, n3):
                v1 = t[i - 1, j - 1, k - 1] + score[mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight \
                     + score[mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight \
                     + score[mapping[seq2[j - 1]], mapping[seq3[k - 1]]]
                v2 = t[i, j - 1, k - 1] + 2 * weighted_gap + score[mapping[seq2[j - 1]], mapping[seq3[k - 1]]]
                v3 = t[i - 1, j, k - 1] + weighted_gap + gap + score[
                    mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight
                v4 = t[i - 1, j - 1, k] + score[
                    mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight + weighted_gap + gap
                v5 = t[i, j, k - 1] + weighted_gap + gap
                v6 = t[i, j - 1, k] + weighted_gap + gap
                v7 = t[i - 1, j, k] + 2 * weighted_gap
                t[i, j, k] = min(v1, v2, v3, v4, v5, v6, v7)
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
                v == t[i - 1, j, k - 1] + weighted_gap + gap + score[mapping[seq1[i - 1]], mapping[seq3[k - 1]]] * weight:
            a1.appendleft(seq1[i - 1])
            a2.appendleft('-')
            a3.appendleft(seq3[k - 1])
            i, k = i - 1, k - 1
        elif i > 0 and j > 0 and \
                v == t[i - 1, j - 1, k] + score[mapping[seq1[i - 1]], mapping[seq2[j - 1]]] * weight + weighted_gap + gap:
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


def sp_score_clique(seqs, clique, k, l):
    """return the sp score of a clique"""
    if l == 2:
        return dynamic_table_2D(seqs[clique[0]], seqs[clique[1]], k - (l - 1))[-1, -1]
    if l == 3:
        return dynamic_table_3D(seqs[clique[0]], seqs[clique[1]], seqs[clique[2]], k - (l - 1))[-1, -1, -1]


def alignment_clique(seqs, clique, k, l):
    """return the optimal alignment of a clique"""
    if l == 2:
        return pairwise_alignment(seqs[clique[0]], seqs[clique[1]], k-(l-1))
    if l == 3:
        return three_exact_alignment(seqs[clique[0]], seqs[clique[1]], seqs[clique[2]], k - (l - 1))


def l_star_align(seqs, l_star, k, l):
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


if __name__ == "__main__":
    names, seqs = parse_fasta("test_seqs/testdata_short.txt")
    a1, a2, a3 = three_exact_alignment(*seqs)
    print(a1, a2, a3)
