import unittest
from decomp import Decomposition
from collections import defaultdict
from math import log, ceil
import ot
from constants import *
from distribution import generate_points
import math
import utils

# TESTS RELATED TO THE CLUSTER DECOMPOSITION
class TestClustering(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.p = P

        # GENERATE THE POINT SET AND DECOMPOSITION
        cls.A, cls.B, cls.masses_A, cls.masses_B = generate_points(N,DIM,DISTRIBUTION)
        distance_function = utils.dist
        decomp = Decomposition(cls.A, cls.B, distance_function)
        cls.clusters = decomp.clusters
        cls.compute_cluster_dist()

    # TEST THAT EVERY PAIR OF POINTS IS REPRESENTED IN SOME CLUSTER
    def test_every_pair_in_some_cluster(self):
        edges = defaultdict(int)
        for center, cluster in self.clusters.items():
            for a in cluster.A:
                for b in cluster.B.keys():
                    edges[(a,b)] = 1

        missing = 0
        for a in self.A:
            for b in self.B:
                if edges[(a,b)] == 0:
                    missing += 1
        self.assertEqual(missing, 0)

    @classmethod
    def compute_cluster_dist(cls):
        # INITIALIZE DICTIONARY FOR CLUSTER DISTANCE
        cls.distC = dict()
        for a in cls.A:
            for b in cls.B:
                cls.distC[(a,b)] = DIAMETER

        # BRUTE FORCE COMPUTE CLUSTER DISTANCE
        for center, cluster in cls.clusters.items():
            for a in cluster.A:
                for b in cluster.B:
                    #if cluster.distC(a,b) < utils.dist(a,b):
                        #print(":(")
                    cls.distC[(a,b)] = min(cluster.distC(a,b), cls.distC[(a,b)])

    def test_clusters_individually_dominate(self):
        oops = 0
        for center, cluster in self.clusters.items():
            for a in cluster.A:
                for b in cluster.B:
                    if cluster.distC(a,b) < utils.dist(a,b):
                        #print("distC: ", cluster.distC(a,b), "dist: ", utils.dist(a,b)) 
                        oops += 1
        self.assertEqual(oops, 0)


    # TEST THAT CLUSTER DISTANCE IS DOMINATING AND APPROXIMATES THE GIVEN METRIC CORRECTLY
    def test_cluster_dist_dict_size(self):
        # SANITY CHECK THAT THE DICTIONARY IS THE CORRECT SIZE
        assert(len(self.distC) == len(self.A) * len(self.B))

    def test_cluster_dist_assigned_to_each_pair(self):
        # CHECK THAT EVERY PAIR IS ASSIGNED A DISTANCE (ie. not 2)
        for a in self.A:
            for b in self.B:
                #if self.distC[(a,b)] == DIAMETER:
                    #print("dist(a,b) = ", utils.dist(a,b), "<? ", DIAMETER)
                self.assertTrue(self.distC[(a,b)] < DIAMETER or DIAMETER <= 2*K*BASE*utils.dist(a,b))

    def test_cluster_dist_is_dominating(self):
        # CHECK CLUSTER DISTANCE IS DOMINATING
        dominates = True
        for a in self.A:
            for b in self.B:
                if self.distC[(a,b)] < utils.dist(a,b):
                    print("distC(a,b) = ", self.distC[(a,b)])
                    print("dist(a,b) = ", utils.dist(a,b))
                    dominates = False
        self.assertTrue(dominates)
                
    def test_approx_factor_of_cluster_dist(self):
        # CHECK APPROXIMATION FACTOR OF CLUSTER DISTANCE
        max_approx = 1
        for a in self.A:
            for b in self.B:
                max_approx = max(max_approx, self.distC[(a,b)]/utils.dist(a,b))
                assert(self.distC[(a,b)] <= 4*BASE*utils.dist(a,b))
        print("\nmax approximation: distC(a,b) <", max_approx, "dist(a,b)")

    # TEST THAT WE DIDN'T ADD ANY CLUSTERS WITH EMPTY A- or B-SETS
    def test_for_empty_A_or_B_clusters(self):
        for center, cluster in self.clusters.items():
            self.assertNotEqual(len(cluster.A), 0)
            self.assertNotEqual(len(cluster.B), 0)

    def test_points_are_in_correct_levels(self):
        for center, cluster in self.clusters.items():
            for a in cluster.A:
                if utils.dist(a,center) == 0:
                    self.assertEqual(cluster.A[a], None)
                else:
                    level = ceil(log(utils.dist(a,center))/log(BASE))
                    self.assertEqual(cluster.A[a], level)
            for b in cluster.B:
                if utils.dist(b,center) == 0:
                    self.assertEqual(cluster.B[b], None)
                else:
                    level = ceil(log(utils.dist(b,center))/log(BASE))
                    self.assertEqual(cluster.B[b], level)

    def test_matching_dist_vs_distC(self):
        dist_matrix = [ [ utils.dist(a,b) for a in self.A] for b in self.B ] #compute_euclidean_distances(self.A, self.B, self.p)
        distC_matrix = [ [ self.distC[(a,b)] for a in self.A ] for b in self.B ]

        cluster_cost = ot.emd2(self.masses_A, self.masses_B, distC_matrix)
        real_cost = ot.emd2(self.masses_A, self.masses_B, dist_matrix)
        ratio = cluster_cost / real_cost
        print("\ndistC matching cost: ", cluster_cost)
        print("dist matching cost: ", real_cost)
        print("ratio: ", ratio)
        self.assertTrue(ratio <= 4*BASE)

if __name__ == "__main__":
    unittest.main()
