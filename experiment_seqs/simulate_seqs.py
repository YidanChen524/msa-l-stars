"""
Simulating fasta files with different number of sequences and length for tests
"""

import random
import os
from helpers import parse_fasta

DNA = 'ACGT'
LINE_WIDTH = 60


def simulate_random_string(m):
    """Simulate a DNA sequence of length m randomly"""
    nucleotides = [random.choice(DNA) for _ in range(m)]
    lines = []
    for i in range(0, m, LINE_WIDTH):
        lines.append(''.join(nucleotides[i:(i + LINE_WIDTH)]))
    return '\n'.join(lines)


def simulate_related_string(template, proportion):
    """Simulate a DNA sequence given a template sequence and maximum edit proportion"""
    m = len(template)
    edits = int(m*proportion)
    nucleotides = list(template)
    for _ in range(edits):
        pos = random.randrange(m)
        nucleotides[pos] = random.choice(list(DNA) + [""])
    lines = []
    for i in range(0, m, LINE_WIDTH):
        lines.append(''.join(nucleotides[i:(i + LINE_WIDTH)]))
    return '\n'.join(lines)


def simulate_random_sequences(n, m):
    """Simulate n sequences of length m."""
    with open(f"{os.path.dirname(os.path.abspath(__file__))}/random_{n}_{m}.fa", "w") as f:
        for i in range(n):
            f.write(f">seq{i}\n")
            f.write(simulate_random_string(m))
            f.write("\n\n")


def simulate_related_sequences(n, m):
    """Simulate n sequences of length m."""
    with open(f"{os.path.dirname(os.path.abspath(__file__))}/related_{n}_{m}.fa", "w") as f:
        template_seq = simulate_random_string(m)
        f.write(f">seq0\n{template_seq}\n\n")
        for i in range(1, n):
            f.write(f">seq{i}\n")
            f.write(simulate_related_string(template_seq, 30))
            f.write("\n\n")


def convert_alignment_to_sequences(file, n, m):
    names, alignment = parse_fasta(file)
    seqs = ["".join([c for c in alm if c != '-']) for alm in alignment]
    with open(f"{os.path.dirname(os.path.abspath(__file__))}/simulated_{n}_{m}.fa", "w") as f:
        for (name, seq) in zip(names, seqs):
            f.write(f">{name}\n")
            f.write(f"{seq}\n")


if __name__ == '__main__':
    for n in (3, 4, 5, 7, 9, 10, 11, 13):
        simulate_random_sequences(n, 10)
