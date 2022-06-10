"""
optimized l-stars algorithm
"""
import sys
from helpers import *
from itertools import combinations


def sp_score_for_all_cliques(seqs, k, l):
    """return a dictionary containing sp score for all possible cliques"""
    scores = {}
    for c in range(k):
        for comb in combinations([i for i in range(k) if i != c], l - 1):
            clique = (c,) + comb
            scores[clique] = sp_score_clique(seqs, clique, k, l)
    return scores


def find_next_cliques(nodes, l, c):
    """yield all available next cliques given a list of nodes which are not added to an l-star"""
    for comb in combinations(nodes[1:], l - 2):
        clique = (c, nodes[0]) + comb
        yield clique


def find_current_collection(nodes, step, l):
    """yield all collections of vertices at the current step"""
    for comb in combinations(nodes[step:], step * (l - 2)):
        yield (*nodes[:step],) + comb


def find_optimal_l_star(seqs, k, l):
    """return the l-star with optimal sp score"""
    # precalculate clique scores
    scores = sp_score_for_all_cliques(seqs, k, l)

    # find optimal l-star for each choice of center
    opt_star, opt_score = None, sys.maxsize
    for c in range(k):
        # pre-fill the dynamic table
        vertices = tuple(v for v in range(k) if v != c)
        dp_table = {}
        for collection in find_current_collection(vertices, 1, l):
            dp_table[collection] = (scores[(c,) + collection], ())

        # fill out the dynamic table step by step
        for step in range(1, (k - 1) // (l - 1)):
            for collection in find_current_collection(vertices, step, l):
                for next_clique in find_next_cliques([v for v in vertices if v not in collection], l, c):
                    new_collection = tuple(sorted(collection + next_clique[1:]))
                    new_score = dp_table[collection][0] + scores[next_clique]
                    if new_collection not in dp_table or dp_table[new_collection][0] > new_score:
                        dp_table[new_collection] = (new_score, collection)

        # backtracking to find an optimal star for center c
        final_score = dp_table[vertices][0]
        if final_score < opt_score:
            l_star = []
            while vertices:
                prev_collection = dp_table[vertices][1]
                l_star.append((c,) + tuple(v for v in vertices if v not in prev_collection))
                vertices = prev_collection
            opt_star, opt_score = l_star, final_score

    return opt_star, opt_score


if __name__ == "__main__":
    # read in sequences and store in Seqs class
    names, seqs = parse_fasta("test_seqs/examples/testdata_7_seqs.txt")
    k, l = len(seqs), 4
    # find the l-star with optimal sp score
    optimal_l_star, optimal_score = find_optimal_l_star(seqs, k, l)
    print(f"Optimal l-star: {optimal_l_star}")
    print(f"Optimal SP-score: {optimal_score}")
    # compute an optimal alignment for the optimal l-star
    alignment = align_l_star(seqs, optimal_l_star, k, l)
    print("Optimal Alignment:")
    print(*alignment, sep="\n")
