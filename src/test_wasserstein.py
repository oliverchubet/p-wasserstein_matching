import unittest
from distribution import generate_points
import utils
from random import randint, seed
from wasserstein import Wasserstein
from point import Point
from math import ceil
from a_points import A, B
from constants import *
import ot

class TestMinimal(unittest.TestCase):

    def test_wasserstein(self):
        for x in [29]: #range(20):
            print("Seed =", x)
            seed(x)
            n = 100
            self.base = 1.01
            self.delta = 0.001
            self.p = 2
            self.A, self.B, masses_A, masses_B = generate_points(n,self.p,"Uniform")
            self.distance_function = utils.dist
            self.wasserstein = Wasserstein(self.A, self.B, self.distance_function, self.p, delta=self.delta, base=self.base)
            if DEBUG: self.wasserstein.print_proxy_dist_matrix()
            print("n =", n, "eps =", self.base-1, "delta =", self.delta, "p =", self.p)
            self.wasserstein.compute_pWasserstein()

            dist_matrix = [ [ pow(utils.dist(a,b),self.p) for a in self.A] for b in self.B ] 
            cluster_dist_matrix = [[ pow(self.wasserstein.distC[(a,b)],self.p) for a in self.A] for b in self.B]
            proxy_dist_matrix = [[ self.wasserstein.proxyDistC[(a,b)] for a in self.A] for b in self.B]
            
            real_cost = pow(ot.emd2(masses_A, masses_B, dist_matrix), 1.0/self.p)
            our_cost = self.wasserstein.min_cost

            our_cluster_cost = self.wasserstein.cost_using_clustering
            emd_cluster_cost = pow(ot.emd2(masses_A, masses_B, cluster_dist_matrix), 1.0/self.p)

            ratio = our_cost / real_cost

            our_cluster_ratio = our_cluster_cost / real_cost
            emd_cluster_ratio = emd_cluster_cost / real_cost

            our_proxy_cost = self.wasserstein.cost_using_proxy
            emd_proxy_cost = pow(ot.emd2(masses_A, masses_B, proxy_dist_matrix), 1.0/self.p)


            print("n =", n, "delta =", self.delta, "eps =", self.base - 1)
            print()

            print("our matching cost: ", our_cost)
            print("\"real\" matching cost: ", real_cost)
            print("ratio: ", ratio)

            print()
            print("our cluster matching cost: ", our_cluster_cost, "\t(cost using our algorithm and cluster distances)")
            print("\"real\" cluster matching cost: ", emd_cluster_cost, "\t(cost using emd and cluster distances)")


            print()
            print("our cluster ratio: ", our_cluster_ratio, "\t(our cluster matching cost / real cost)")
            print("\"real\" cluster ratio: ", emd_cluster_ratio, "\t(\"real\" cluster matching cost / real cost)")


            print()
            print("Note: ratio should be less than our cluster ratio, or something is wrong")
            print("Note: our cluster ratio should match the \"real\" cluster ratio for small enough delta")

            print("our proxy matching cost:", our_proxy_cost)
            print("\"real\" proxy cost:", emd_proxy_cost)

            #print([[ceil(pow(x,2)/self.delta) for x in row] for row in cluster_dist_matrix])

