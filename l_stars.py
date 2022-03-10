"""
l-stars with dynamic programming
"""
import sys
from helpers import *
from itertools import combinations


def generate_all_l_stars(k, l):
    """given the length of the sequences k and size of each clique l, return all possible l-stars"""
    l_stars = []

    def _next_clique(center, nodes):
        for comb in combinations(nodes[1:], l-2):
            clique = (center, nodes[0],) + comb
            remain = [i for i in nodes if i not in clique]
            yield clique, remain

    def _backtrack(nodes, center, cliques):
        if not nodes:
            l_stars.append(cliques.copy())
        else:
            for clique, remain in _next_clique(center, nodes):
                cliques.append(clique)
                _backtrack(remain, center, cliques)
                cliques.pop()

    for c in range(k):
        _backtrack([i for i in range(k) if i != c], c, [])

    return l_stars


def sp_score_for_all_cliques(seqs, k, l):
    """return a dictionary containing sp score for all possible cliques"""
    scores = {}
    for c in range(k):
        for comb in combinations([i for i in range(k) if i != c], l-1):
            clique = (c,) + comb
            scores[clique] = sp_score_clique(seqs, clique, k, l)
    return scores


def find_best_l_star(seqs, k, l):
    """return the l-star with optimal sp score"""
    l_stars, scores = generate_all_l_stars(k, l), sp_score_for_all_cliques(seqs, k, l)
    opt_score, opt_l_star = sys.maxsize, None
    for l_star in l_stars:
        s = 0
        for clique in l_star:
            s += scores[clique]
        if s < opt_score:
            opt_score, opt_l_star = s, l_star
    return opt_l_star, opt_score


if __name__ == "__main__":
    # read in sequences and store in Seqs class
    names, seqs = parse_fasta("test_seqs/testdata_7_seqs.txt")
    k, l = len(seqs), 2
    # find the l-star with optimal sp score
    optimal_l_star, optimal_score = find_best_l_star(seqs, k, l)
    print(f"Optimal l-star: {optimal_l_star}")
    print(f"Optimal SP-score: {optimal_score}")
    # compute an optimal alignment for the optimal l-star
    alignment = l_star_align(seqs, optimal_l_star, k, l)
    print("Optimal Alignment:")
    print(*alignment, sep="\n")

