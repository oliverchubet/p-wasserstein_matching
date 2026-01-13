# Computing p-Wasserstein Distance Using Min-Cost Matching

Inputs:
    - A, B: two n-point sets in some metric space X
    - cost: a metric on A x B
    - p: integer (specifying the Wp-metric)
Parameters:
    - delta: additive approximation
    - eps: multiplicative approximation
Output:
    - M: a matching where the weight is a (1+eps)-multiplicative, delta-additive approximation of Wp(A,B)

Algorithm:
1. Initialization
    - Initialize M to an empty matching
    - Initialize dual weights y(x) = 0 for each point in A and B
    - Compute a BCP data structure on A x B
2. Phase loop
    - Djikstra: Compute the shortest path using BCP
    - Dual weights: Update the dual weights using the length of the shortest path
    - Partial DFS: Find vertex disjoint set of augmenting paths
    - Augmentation: Update M using the paths found
3. Compute and return the cost
