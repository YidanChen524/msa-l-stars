"""
randomized l-stars algorithm
"""
import sys
import math
from random import sample
from helpers import *


def randomized_l_star(k, l, center):
    """generate a randomized l-star"""
    random_sample = sample([i for i in range(k) if i != center], k-1)
    l_star = []
    for i in range(0, k-1, l-1):
        l_star.append((center,) + tuple(random_sample[i:i+l-1]))
    return l_star


def find_optimal_randomized_l_star(seqs, k, l, epsilon):
    """find the optimal l-star returned by the randomized algorithm"""
    opt_score, opt_star = sys.maxsize, None
    for c in range(k):
        for _ in range(int(2 * math.log(k / epsilon, 2))):
            l_star = randomized_l_star(k, l, c)
            tmp_score = sum([sp_score_clique(seqs, clique, k, l) for clique in l_star])
            if tmp_score < opt_score:
                opt_score = tmp_score
                opt_star = l_star
    return opt_star, opt_score
