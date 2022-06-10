"""
Experiments
"""

import time
from helpers import sp_score, parse_fasta, align_l_star, align_2l_star
from optimized_l_stars import find_optimal_l_star
from paired_l_stars import find_optimal_star
from randomized_l_stars import find_optimal_randomized_l_star
from test_seqs.simulate_seqs import simulate_random_sequences
import csv


def evaluation(func):
    """print the runtime and sp score of the given algorithm"""

    def wrapper_time_score(*args, **kwargs):
        start_time = time.perf_counter()
        alignment = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        score = sp_score(alignment)
        print(f"file {kwargs['file']} using {func.__name__} with l={kwargs['l']} time: {run_time:.4f}, score: {score}")
        return run_time, score

    return wrapper_time_score


@evaluation
def optimized_l_star(file, k, l):
    _, seqs = parse_fasta(file)
    opt_star, _ = find_optimal_l_star(seqs, k, l)
    return align_l_star(seqs, opt_star, k, l)


@evaluation
def paired_l_stars(file, k, l):
    _, seqs = parse_fasta(file)
    opt_star, _ = find_optimal_star(seqs, k, l)
    if len(opt_star[0]) == l:
        return align_l_star(seqs, opt_star, k, l)
    else:
        return align_2l_star(seqs, opt_star, k, l)


@evaluation
def randomized_l_stars(file, k, l, eps):
    _, seqs = parse_fasta(file)
    opt_star, _ = find_optimal_randomized_l_star(seqs, k, l, eps)
    return align_l_star(seqs, opt_star, k, l)


def test_approximation_ratios():
    """tests the approximation ratios"""
    pass


def test_optimized_l_stars():
    """
    l = 2, 3
    k = 5, 7, 9, 13
    """
    pass


def test_paired_l_stars():
    """
    l = 2, 3
    k = 5, 7, 9, 13
    """
    pass


def test_randomized_l_stars():
    """
    l = 2
    k = 4, 7, 10, 13
    """
    pass


if __name__ == "__main__":
    # l-stars-method
    # for k in (7, 13):
    #     for m in (10,):
    #         for l in (2, 3, 4):
    #             l_stars_method(file=f"test_seqs/experiments/random_{k}_{m}.fa", k=k, l=l)
    # (2l-1)-stars
    # for k in (9,):
    #     for m in (10,):
    #         for _ in range(3):
    #             simulate_random_sequences(k, m)
    #             for l in (2, 3):
    #                 l_stars_small_balanced_set(file=f"test_seqs/experiments/random_{k}_{m}.fa", k=k, l=l)
    # randomized algorithm
    # for k in (7, 13, 19):
    #     for m in (10,):
    #         for l in (2, 3, 4):
    #             randomized_l_stars(file=f"test_seqs/experiments/random_{k}_{m}.fa", k=k, l=l, eps=0.01)
    # k=13 for all algo
    simulate_random_sequences(13, 10)
    file = f"test_seqs/experiments/random_13_10.fa"
    # l_stars_method(file=file, k=9, l=2)
    # l_stars_method(file=file, k=9, l=3)
    paired_l_stars(file=file, k=13, l=2)
    paired_l_stars(file=file, k=13, l=3)
    # r = 100
    # t1, s1 = [None] * r, [None] * r
    # t2, s2 = [None] * r, [None] * r
    # t3, s3 = [None] * r, [None] * r
    # t4, s4 = [None] * r, [None] * r
    # t5, s5 = [None] * r, [None] * r
    # for i in range(r):
    #     t1[i], s1[i] = randomized_l_stars(file=file, k=13, l=3, eps=0.1)
    #     t2[i], s2[i] = randomized_l_stars(file=file, k=13, l=3, eps=0.3)
    #     t3[i], s3[i] = randomized_l_stars(file=file, k=13, l=3, eps=0.5)
    #     t4[i], s4[i] = randomized_l_stars(file=file, k=13, l=3, eps=0.7)
    #     t5[i], s5[i] = randomized_l_stars(file=file, k=13, l=3, eps=0.9)
    # print(f"eps=0.1: time = {sum(t1)/r}, score = {sum(s1)/r}, 90th percentile = {sorted(s1)[int(r*0.9)-1]}")
    # print(f"eps=0.3: time = {sum(t2)/r}, score = {sum(s2)/r}, 70th percentile = {sorted(s2)[int(r*0.7)-1]}")
    # print(f"eps=0.5: time = {sum(t3)/r}, score = {sum(s3)/r}, 50th percentile = {sorted(s3)[int(r*0.5)-1]}")
    # print(f"eps=0.7: time = {sum(t4)/r}, score = {sum(s4)/r}, 30th percentile = {sorted(s4)[int(r*0.3)-1]}")
    # print(f"eps=0.9: time = {sum(t5)/r}, score = {sum(s5)/r}, 10th percentile = {sorted(s5)[int(r*0.1)-1]}")

    # with open("test_results/randomized_eps.csv", "w") as wf:
    #     writer = csv.writer(wf)
    #     writer.writerow(["eps", "time", "score"])
    #     [writer.writerow(row) for row in zip([0.1] * r, t1, s1)]
    #     [writer.writerow(row) for row in zip([0.3] * r, t2, s2)]
    #     [writer.writerow(row) for row in zip([0.5] * r, t3, s3)]
    #     [writer.writerow(row) for row in zip([0.7] * r, t4, s4)]
    #     [writer.writerow(row) for row in zip([0.9] * r, t5, s5)]

