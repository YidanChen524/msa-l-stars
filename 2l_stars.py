"""
(2l-1)-stars algorithms
"""
from helpers import *


def generate_l_star(k, l, center):
    """given k, l and the center, generate an arbitrary l-star"""
    l_star = []
    clique_nodes = [i for i in range(k) if i != center]
    for i in range(0, k-1, l-1):
        l_star.append((center,) + tuple(clique_nodes[i:i+l-1]))
    return l_star


def graph(seqs, l_star, k , l):
    """given sequences and an l-star, return the corresponding graph"""
    n = len(l_star)
    g = np.zeros([n, n])
    for i in range(n):
        for j in range(i+1, n):
            g[i, j] = g[j, i] = sp_score_clique_2l_star(seqs, l_star[i]+l_star[j][1:], k, l)
    return g


def optimal_2l_star(g):
    """given a graph g, return the optimal (2l-1)-star by solving the matching problem"""
    pass


if __name__ == "__main__":
    # read in sequences and store in Seqs class
    names, seqs = parse_fasta("test_seqs/testdata_7_seqs.txt")
    k, l = len(seqs), 2
    l_star = generate_l_star(k, l, 0)
    print(l_star)
    g = graph(seqs, l_star, k, l)
    print(g)
