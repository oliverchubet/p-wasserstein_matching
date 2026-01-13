from constants import *
from distribution import generate_points
from utils import dist
from ai_decomp import Decomposition
from itertools import chain
from tree import Tree
import copy
from math import ceil, log
import ot
import numpy as np
from scipy.spatial.distance import cdist

# Assume constants, generate_points, etc. are imported as before

#class WassersteinOptimized:
class Wasserstein:
    def __init__(self, A, B, dist_func, p=1, delta=0.01, base=1.01, dim=2):
        self.A = A
        self.B = B
        self.n = len(A)
        self.p = p
        
        # MAP OBJECTS TO INDICES for array access
        self.a_to_idx = {a: i for i, a in enumerate(A)}
        self.b_to_idx = {b: i for i, b in enumerate(B)}
        
        # Pre-compute raw coordinates for vectorization
        # Assuming your objects a/b have .coords or similar, or are tuples
        # If A is a list of objects with x,y properties:
        #self.A_coords = np.array([[a.x, a.y] for a in A]) 
        #self.B_coords = np.array([[b.x, b.y] for b in B])
        self.A_coords = np.array([a.coords for a in A]) 
        self.B_coords = np.array([b.coords for b in B])
        
        self.diam = 2 * pow(2, dim)
        
        # OPTIMIZATION 1: Use Numpy Arrays instead of Dicts
        # Initialize with diameter
        self.distC = np.full((self.n, self.n), self.diam, dtype=np.float64)
        self.proxyDistC = np.full((self.n, self.n), np.ceil(self.diam/delta), dtype=np.float64)
        
        # Keep your decomposition logic (assuming it accepts the objects)
        self.decomp = Decomposition(A, B, dist_func, p, delta, k=2, base=base, dim=dim)
        
        # Set management
        self.freeB = set(self.B)
        self.matched = 0
    
    def compute_cluster_distance(self):
        """
        Vectorized update of cluster distances.
        This replaces the O(N^2) dictionary loop.
        """
        # If Decomposition allows iterating clusters, we update the matrix blocks
        # This part depends on how 'cluster' is structured. 
        # Ideally, 'cluster.A' returns a list of indices, not objects.
        
        for center, cluster in self.decomp.clusters.items():
            # Get indices of points in this cluster
            # This is a hypothetical optimization: 
            # You should implement .indices properties in your cluster objects
            a_indices = [self.a_to_idx[a] for a in cluster.A]
            b_indices = [self.b_to_idx[b] for b in cluster.B]
            
            if not a_indices or not b_indices:
                continue
                
            # Create a meshgrid for the indices to update the block
            # Note: This logic mimics your loops but does it in blocks if possible.
            # If your clusters are single points, this doesn't help much. 
            # If they are groups, this is huge.
            
            # For now, if we must loop because cluster.distC is a function:
            for a in cluster.A:
                i = self.a_to_idx[a]
                for b in cluster.B:
                    j = self.b_to_idx[b]
                    
                    # Direct array access is faster than dict hashing
                    d_val = cluster.distC(a, b)
                    if d_val < self.distC[i, j]:
                        self.distC[i, j] = d_val
                        
                    p_val = cluster.proxyDistC(a, b)
                    if p_val < self.proxyDistC[i, j]:
                        self.proxyDistC[i, j] = p_val

    def compute_cost(self):
        """
        Vectorized cost computation
        """
        # Gather matches
        matches = [(a, a.match) for a in self.A if a.match]
        if not matches:
            return 0
            
        a_objs, b_matched = zip(*matches)
        
        # Convert to indices
        idx_a = [self.a_to_idx[a] for a in a_objs]
        idx_b = [self.b_to_idx[b] for b in b_matched]
        
        # Vectorized lookup from the precomputed matrix
        # This replaces the sum([pow(...)]) loop
        costs = self.distC[idx_a, idx_b]
        total_cost = np.sum(np.power(costs, self.p))
        self.cost_using_clustering = pow(total_cost, 1.0/self.p)
        
        # Real Euclidean dist calculation
        # Assuming 2D points
        real_dists = np.linalg.norm(self.A_coords[idx_a] - self.B_coords[idx_b], axis=1)
        self.min_cost = pow(np.sum(np.power(real_dists, self.p)), 1.0/self.p)

# ---------------------------------------------------------
# OPTIMIZED MAIN EXECUTION
# ---------------------------------------------------------
if __name__ == "__main__":
    n = 100
    p = 2
    
    # 1. Generate Points
    A, B, masses_A, masses_B = generate_points(n, p, "Uniform")
    
    # Extract coordinates into numpy arrays immediately for POT and Verification
    # (Assuming A and B are objects with x, y attributes or similar)
    # If A/B are just tuples, use: A_coords = np.array(A)
    #A_coords = np.array([[a.x, a.y] for a in A]) 
    #B_coords = np.array([[b.x, b.y] for b in B])
    A_coords = np.array([a.coords for a in A]) 
    B_coords = np.array([b.coords for b in B])
    
    # 2. Compute Ground Truth Matrix (Vectorized)
    # This replaces: [ [ pow(dist(a,b),2) ... ] ... ]
    # cdist computes the distance matrix roughly 100x faster than list comp
    print("Computing Ground Truth Matrix...")
    dist_matrix = cdist(A_coords, B_coords, metric='euclidean')
    dist_matrix = np.power(dist_matrix, p) # Apply power p
    
    # 3. Run Your Algorithm
    print("Running Wasserstein...")
    wasserstein = Wasserstein(A, B, dist, p, delta=0.01, base=1.01)
    #wasserstein.compute_pWasserstein()
    wasserstein.compute_cost()
    
    # 4. Vectorized Cost Comparison
    # We map the dict back to a matrix for `ot.emd2`
    # (The optimized class above would already have this as an array)
    cluster_dist_matrix = np.zeros((n, n))
    
    # Fast fill of the cluster matrix from your existing dictionary structure
    # (Only needed if you keep the dictionary implementation)
    # Since we can't vectorize dict access easily, we just iterate once efficiently
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            cluster_dist_matrix[i, j] = pow(wasserstein.distC.get((a,b), 0), 2)

    # 5. Compute POT metrics
    real_cost = pow(ot.emd2(masses_A, masses_B, dist_matrix), 1.0/p)
    cluster_emd_cost = pow(ot.emd2(masses_A, masses_B, cluster_dist_matrix), 1.0/p)
    
    print(f"Our matching cost: {wasserstein.min_cost}")
    print(f"Real matching cost: {real_cost}")
