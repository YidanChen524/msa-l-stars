"""
Experiments
"""

import time
from helpers import sp_score, parse_fasta, align_l_star, align_2l_star
from optimized_l_stars import find_optimal_l_star
from paired_l_stars import find_optimal_star
from randomized_l_stars import find_optimal_randomized_l_star
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
def optimized_l_stars(file, k, l):
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


def test_optimized_l_stars(round):
    """
    l = 2, k = 3, 5, 7, 9, 11, 13
    l = 3, k = 3, 5, 7, 9, 11, 13
    l = 4, k = 4, 7, 10, 13
    """
    with open("experiment_results/optimized_l_stars.csv", "a") as wf:
        writer = csv.writer(wf)
        for k in (3, 5, 7, 9, 11, 13):
            for l in (2, 3):
                rt, sc = optimized_l_stars(file=f"experiment_seqs/round_{round}/random_{k}_10.fa", k=k, l=l)
                writer.writerow([k, l, rt, sc, round])
        for k in (4, 7, 10, 13):
            for l in (4,):
                rt, sc = optimized_l_stars(file=f"experiment_seqs/round_{round}/random_{k}_10.fa", k=k, l=l)
                writer.writerow([k, l, rt, sc, round])


def test_paired_l_stars(round):
    """
    l = 2, k = 3, 5, 7, 9, 11, 13
    l = 3, k = 5, 9, 13
    """
    with open("experiment_results/paired_l_stars.csv", "a") as wf:
        writer = csv.writer(wf)
        for k in (3, 5, 7, 9, 11, 13):
            for l in (2, ):
                rt, sc = paired_l_stars(file=f"experiment_seqs/round_{round}/random_{k}_10.fa", k=k, l=l)
                writer.writerow([k, l, rt, sc, round])
        for k in (5, 9, 13):
            for l in (3,):
                rt, sc = paired_l_stars(file=f"experiment_seqs/round_{round}/random_{k}_10.fa", k=k, l=l)
                writer.writerow([k, l, rt, sc, round])


def test_randomized_l_stars(round, eps=0.1):
    """
    l = 2, k = 3, 5, 7, 9, 11, 13
    l = 3, k = 3, 5, 7, 9, 11, 13
    l = 4, k = 4, 7, 10, 13
    """
    with open("experiment_results/randomized_l_stars.csv", "a") as wf:
        writer = csv.writer(wf)
        for k in (3, 5, 7, 9, 11, 13):
            for l in (2, 3):
                rt, sc = randomized_l_stars(file=f"experiment_seqs/round_{round}/random_{k}_10.fa", k=k, l=l, eps=eps)
                writer.writerow([k, l, rt, sc, round, eps])
        for k in (4, 7, 10, 13):
            for l in (4,):
                rt, sc = randomized_l_stars(file=f"experiment_seqs/round_{round}/random_{k}_10.fa", k=k, l=l, eps=eps)
                writer.writerow([k, l, rt, sc, round, eps])


if __name__ == "__main__":
    for r in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10):
        test_optimized_l_stars(r)
        test_paired_l_stars(r)
        for eps in (0.1, 0.3, 0.6, 0.9):
            test_randomized_l_stars(r, eps)
