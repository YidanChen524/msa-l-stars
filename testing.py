"""
Experiments:
    - test if generalized l-stars method has better performance (l = 2, 3, 4, 5)
    - test if using a smaller balanced set is faster and its performance
    - test if randomized algorithm meets requirements and its performance
"""

import time
from helpers import sp_score, parse_fasta, align_l_star, align_2l_star
from l_stars import find_optimal_l_star
from l_stars_small_balanced_set import find_optimal_2l_star
from randomized_l_stars import find_optimal_randomized_l_star


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
def l_stars_method(file, k, l):
    _, seqs = parse_fasta(file)
    opt_star, _ = find_optimal_l_star(seqs, k, l)
    return align_l_star(seqs, opt_star, k, l)


@evaluation
def l_stars_small_balanced_set(file, k, l):
    _, seqs = parse_fasta(file)
    opt_star, _ = find_optimal_2l_star(seqs, k, l)
    return align_2l_star(seqs, opt_star, k, l)


@evaluation
def randomized_l_stars(file, k, l, eps):
    _, seqs = parse_fasta(file)
    opt_star, _ = find_optimal_randomized_l_star(seqs, k, l, eps)
    return align_l_star(seqs, opt_star, k, l)


if __name__ == "__main__":
    for k in (9, 13, 17):
        for m in (10,):
            for l in (2, 3):
                l_stars_small_balanced_set(file=f"test_seqs/experiments/related_{k}_{m}.fa", k=k, l=l)

