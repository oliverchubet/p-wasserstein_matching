# README

## Testing

The code can be tested with `python -m unittest`.

## Computing the Cluster Decomposition

There are four main steps to construct the (2-layer) cluster decomposition.
    1. Compute the sample.
    2. For each point in the sample, compute a cluster.
    3. For each input point (all points) compute the distance to the closest point in the sample.
    4. For each non-sample point, compute a cluster.

### Computing the Sample

For an input of size $n$, each point is added to the sample with probability $1/\sqrt{n}$.

### Computing Clusters for Sample Points

The clusters of sample points contain _all_ points.
To construct a cluster, we iterate through all points and put each in a bucket keyed by level of the distance to the center, computed as:
        `level = ceil(log(self.dist(a,center),BASE))`,
where BASE is a user-specified approximation factor.
Note that A-points and B-points each have their own set of buckets, `bucketsA` and `bucketsB`.
The points are passed in two dictionaries `pointsA` and `pointsB` that map each point to its level.
For each point, we also add the center to a reverse dictionary `lookup` mapping points to the clusters that contain it.

### Dictionary of Nearest Neighbor Distances

The nearest neighbor distances from each point to the sample are computed in a brute-force fashion.
We iterating through each point in the sample for each point,and take the smallest distance found.

### Computing Clusters for Non-Sample Points

The clusters of non-sample points only contain those points whose distance to the center is less than the nearest neighbor distance to the sample.
In all other respects, the computation of the clusters for non-sample points is the same as for sample points.

## Dynamic Weighted Bichromatic Closest Pair

We use the decomposition to compute weighted bichromatic closest pair.
    1. Each cluster constructs max-heaps for each non-empty bucket and computes an overall min-slack edge for the cluster.
    2. We add the min-slack edge of each cluster to a min-heap.

### Computing the Min-Slack Edge of a Cluster

Each point adds its weight to the max-heap for its level.
For our purposes we assume that the points of lower levels are also contained in the larger levels, so the max weight of a given level i will be be the max weight among all levels at most i.
Then, to compute the min-slack edge we loop through each non-empty level from smallest to largest and keep track of the max of the max weights seen `max_a` and `max_b` and their corresponding points `a` and `b`, as well as the min-slack and min-slack edge.
This should take at most log(spread) time per cluster.

However, we can ammortize the cost over all clusters.
In expectation each cluster participates in sqrt(n) clusters, so the total number of levels over all clusters is $n^{3/2}$ in expectation.
Within a cluster, there are $O(1)$ heap operations for a (non-empty) level. 
So, computing the min-slack edge for all clusters is $O(n^{3/2}\log(n))$.

### Updating the Weight of a Point

To update the weight of a point or remove a point, we use the reverse dictionary to iterate through the relevant clusters.
Each of those clusters will update their max heap and recompute the min-slack edge for that cluster.

In expectation $\sqrt{n}$ clusters are updated.
The max-heap update is one heap operation.
To recompute the min-slack edge for each cluster is at most log(spread) time.
Each cluster affected also changes priority in the min-heap.

In total, the update time is $\sqrt{n}(\log(n) + \log(\Delta))$ in expectation.

