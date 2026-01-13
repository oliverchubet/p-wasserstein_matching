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
        #for count in range(100):
        # error for n=3, seed(54
        #for x in range(100):
        for x in [54]:
            print("Seed =", x)
            seed(x)
            n = 5000
            self.base = 1.01
            self.delta = 0.01
            self.p = 2
            self.A, self.B, masses_A, masses_B = generate_points(n,self.p,"Uniform")
            self.distance_function = utils.dist
            self.wasserstein = Wasserstein(self.A, self.B, self.distance_function, self.p, delta=self.delta, base=self.base)
            if DEBUG: self.wasserstein.print_proxy_dist_matrix()
            self.wasserstein.compute_pWasserstein()

            #dist_matrix = [ [ pow(utils.dist(a,b),2) for a in self.A] for b in self.B ] 
            dist_matrix = [ [ pow(utils.dist(a,b),2) for a in self.A] for b in self.B ] 
            #cluster_dist_matrix = [[ ceil(pow(self.wasserstein.distC[(a,b)],2)/self.delta)*self.delta for a in self.A] for b in self.B]
            cluster_dist_matrix = [[ pow(self.wasserstein.distC[(a,b)],2) for a in self.A] for b in self.B]
            
            real_cost = pow(ot.emd2(masses_A, masses_B, dist_matrix), 1.0/2)
            cost = self.wasserstein.min_cost
            cluster_cost = self.wasserstein.cost_using_clustering
            cluster_emd_cost = pow(ot.emd2(masses_A, masses_B, cluster_dist_matrix), 1.0/2)
            cluster_ratio = cluster_cost / real_cost
            cluster_emd_ratio = cluster_emd_cost / real_cost
            ratio = cost / real_cost
            print()

            print("our matching cost: ", cost)
            print("\"real\" matching cost: ", real_cost)
            print("ratio: ", ratio)

            print()
            print("our cluster matching cost: ", cluster_cost)
            print("\"real\" matching cost: ", cluster_emd_cost)

            print()
            print("our cluster ratio: ", cluster_ratio)
            print("\"real\" cluster ratio: ", cluster_emd_ratio)


