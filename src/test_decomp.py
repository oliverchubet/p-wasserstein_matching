import unittest
from decomp import Decomposition
from collections import defaultdict
from math import log, ceil
import ot
from constants import *
from distribution import generate_points
import math
import utils

distance_function = utils.dist

class TestClustering(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.p = 2
        cls.delta = 0.01
        cls.k = 2
        cls.n = 10
        cls.dim = 2
        cls.base = 1.01
        cls.distribution = "Normal"
        cls.A, cls.B, cls.masses_A, cls.masses_B = generate_points(cls.n,cls.dim,cls.distribution)
        cls.decomp = Decomposition(cls.A, cls.B, distance_function, cls.p, cls.delta, cls.k, cls.base, cls.dim)
        cls.clusters = cls.decomp.clusters
        cls.decomp.compute_cluster_dist()
        cls.distC = cls.decomp.distC
        cls.diam = cls.decomp.diam

    # TEST THAT EVERY PAIR OF POINTS IS REPRESENTED IN SOME CLUSTER
    def test_every_pair_in_some_cluster(self):
        edges = defaultdict(int)
        for center, cluster in self.clusters.items():
            for a in cluster.A:
                for b in cluster.B:
                    edges[(a,b)] = 1
        missing = 0
        for a in self.A:
            for b in self.B:
                if edges[(a,b)] == 0:
                    missing += 1
        self.assertEqual(missing, 0)

    def test_clusters_individually_dominate(self):
        oops = 0
        for center, cluster in self.clusters.items():
            for a in cluster.A:
                for b in cluster.B:
                    if self.distC[(a,b)] < distance_function(a,b):
                        oops += 1
        self.assertEqual(oops, 0)


    def test_cluster_dist_dict_size(self):
        assert(len(self.distC) == len(self.A) * len(self.B))

    def test_cluster_dist_assigned_to_each_pair(self):
        for a in self.A:
            for b in self.B:
                self.assertTrue(self.distC[(a,b)] < self.diam or self.diam <= 2*self.k*self.base*distance_function(a,b))

    def test_cluster_dist_is_dominating(self):
        dominates = True
        for a in self.A:
            for b in self.B:
                if self.distC[(a,b)] < distance_function(a,b):
                    dominates = False
        self.assertTrue(dominates)
                
    def test_approx_factor_of_cluster_dist(self):
        max_approx = 1
        for a in self.A:
            for b in self.B:
                max_approx = max(max_approx, self.distC[(a,b)]/utils.dist(a,b))
                assert(self.distC[(a,b)] <= 4*self.base*distance_function(a,b))

    def test_for_empty_A_or_B_clusters(self):
        for center, cluster in self.clusters.items():
            self.assertNotEqual(len(cluster.A), 0)
            self.assertNotEqual(len(cluster.B), 0)

    def test_sanity_log_works_as_expected(self):
        self.assertEqual(ceil(log(0.032)/log(1.01)), ceil(log(0.032, 1.01)))

    def test_points_are_in_correct_levels(self):
        for center, cluster in self.clusters.items():
            for a in cluster.A:
                if distance_function(a,center) == 0:
                    self.assertEqual(cluster.A[a], None)
                else:
                    level = ceil(log(distance_function(a,center),self.base))
                    self.assertEqual(cluster.A[a], level)
            for b in cluster.B:
                if distance_function(b,center) == 0:
                    self.assertEqual(cluster.B[b], None)
                else:
                    level = ceil(log(distance_function(b,center),self.base))
                    self.assertEqual(cluster.B[b], level)

    def test_matching_dist_vs_distC(self):
        dist_matrix = [ [ distance_function(a,b) for a in self.A] for b in self.B ]
        distC_matrix = [ [ self.distC[(a,b)] for a in self.A ] for b in self.B ]

        cluster_cost = ot.emd2(self.masses_A, self.masses_B, distC_matrix)
        real_cost = ot.emd2(self.masses_A, self.masses_B, dist_matrix)
        ratio = cluster_cost / real_cost
        #print("emd distC matching cost: ", cluster_cost)
        #print("emd dist matching cost: ", real_cost)
        #print("ratio: ", ratio)
        self.assertTrue(ratio <= 2*self.k*self.base)

if __name__ == "__main__":
    unittest.main()
