"""
(2l-1)-stars algorithms
"""
import sys
from helpers import *
import networkx as nx


def generate_l_star(k, l, center):
    """given k, l and the center, generate an arbitrary l-star"""
    l_star = []
    clique_nodes = [i for i in range(k) if i != center]
    for i in range(0, k-1, l-1):
        l_star.append((center,) + tuple(clique_nodes[i:i+l-1]))
    return l_star


def graph(seqs, l_star, k, l):
    """given sequences and an l-star, return the corresponding graph"""
    n = len(l_star)
    g = np.zeros([n, n])
    for i in range(n):
        for j in range(i+1, n):
            g[i, j] = g[j, i] = sp_score_clique_2l_star(seqs, l_star[i]+l_star[j][1:], k, l)
    return g


def find_optimal_star(seqs, k, l):
    """find the optimal 2l-1 star by iterating through all center strings"""
    opt_score, opt_star = sys.maxsize, None
    for c in range(k):
        # score of the chosen arbitrary l star
        l_star = generate_l_star(k, l, c)
        l_star_score = sum(sp_score_clique(seqs, clique, k, l) for clique in l_star)
        if l_star_score < opt_score:
            opt_score, opt_star = l_star_score, l_star
        # find optimal (2l-1)-star
        g = graph(seqs, l_star, k, l)
        m = nx.min_weight_matching(nx.Graph(g))
        temp_score = sum([g[a, b] for (a, b) in m])
        if temp_score < opt_score:
            opt_score = temp_score
            opt_star = [l_star[a] + l_star[b][1:] for (a, b) in m]
    return opt_star, opt_score
