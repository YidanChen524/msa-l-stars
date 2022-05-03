# Extended Approximation Algorithms
Extended approximation algorithms for multiple sequence alignment

## Introduction
- what is multiple sequence alignment and what are our goal?
- Exact Algorithm.
- Approximation Algorithms.

## L-star methods

### Definitions
- What is it?
- How to compute an alignment based on l-star?
- How and why it can be used to find the optimal alignment with an approximation ratio of 2-l/k

### L-stars with dynamic programming
- l = 2, 3
- Typical Dynamic programming
- Optimization

### (2l-1)-stars
- l = 2, 3
- generate random l-star for each center
- calculate score for each pair
- solve the matching problem
- return the best l-star and best alignment

### random sampling


## Implementation
## Experiments
goals:
- show generalized l-star methods improves quality of alignment when l increases
- 2l-star and randomized methods improve running time while guaranteed the quality of the alignment

## Conclusion

## To Do
- generate sequences with known sp score
- gather data on time and score of each algorithm
- visualization
- when k-1 cannot be divided by l-1

## Questions:
- what kind of biological sequences suitable for l-stars?
- sequences with known sp score, how to prove the approximation ratio is met
- prove the running time in practice
- format of the thesis

# Results


