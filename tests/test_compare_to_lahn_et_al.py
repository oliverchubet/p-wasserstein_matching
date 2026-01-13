import unittest
from constants import *
from distribution import generate_points
import utils
from hungarian import Hungarian
import ot
from math import ceil

class TestApprox(unittest.TestCase):

    def test_min_cost_matching_Uniform(self):
        for i in range(1):
            for n in [10]:
                print("Uniform Test, ",n,"points")
                self.delta = 1/(10*n)
                self.A, self.B, self.masses_A, self.masses_B = generate_points(n,2,"Uniform")
                distance_function = utils.dist #lambda x,y: utils.rounded_dist(x,y,delta=delta)
                self.hungarian = Hungarian(self.A, self.B, distance_function, 2, delta=self.delta)
                self.hungarian.reset_matching()
                self.hungarian.compute_min_cost_matching()
                dist_matrix = [ [ pow(utils.dist(a,b),2) for a in self.A] for b in self.B ] 
                cluster_dist_matrix = [[ ceil(pow(self.hungarian.distC[(a,b)],2)/self.delta)*self.delta for a in self.A] for b in self.B]
                real_cost = pow(ot.emd2(self.masses_A, self.masses_B, dist_matrix), 1.0/2)
                cluster_emd_cost = pow(ot.emd2(self.masses_A, self.masses_B, cluster_dist_matrix), 1.0/2)
                cluster_cost = self.hungarian.cost_using_clustering
                cost = self.hungarian.min_cost
                cluster_ratio = cluster_cost / real_cost
                cluster_emd_ratio = cluster_emd_cost / real_cost
                ratio = cost / real_cost
                print()
                print("our matching cost: ", cost)
                print("cluster matching cost: ", cluster_cost)
                print("cluster emd matching cost: ", cluster_emd_cost)
                print("emd2 matching cost: ", real_cost)
                print("ratio: ", ratio)
                print("cluster ratio: ", cluster_ratio)
                print("cluster emd ratio: ", cluster_emd_ratio)

                with open("1_10n_delta", "a") as file:
                    file.write("\t".join([str(n), "2", "2", str(self.delta), str(cost), str(cluster_cost), str(real_cost), str(ratio),str(cluster_ratio),str(cluster_emd_ratio),str(self.hungarian.edges_seen),str(self.hungarian.phases),"\n"]))



