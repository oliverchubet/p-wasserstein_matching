import unittest
from constants import *
from distribution import generate_points
import utils
from hungarian import Hungarian
import ot

class TestApprox(unittest.TestCase):

    def test_min_cost_matching_Normal(self):
        for i in range(1):
            print("Normal test ", i)
            self.A, self.B, self.masses_A, self.masses_B = generate_points(N,DIM,"Normal_same")
            distance_function = utils.dist
            self.hungarian = Hungarian(self.A, self.B, distance_function, P)
            self.hungarian.reset_matching()
            self.hungarian.compute_min_cost_matching()
            dist_matrix = [ [ pow(utils.dist(a,b),P) for a in self.A] for b in self.B ] 
            real_cost = pow(ot.emd2(self.masses_A, self.masses_B, dist_matrix), 1.0/P)
            cost = self.hungarian.min_cost
            ratio = cost / real_cost
            print()
            print("our matching cost: ", cost)
            print("emd2 matching cost: ", real_cost)
            print("ratio: ", ratio)

    def test_min_cost_matching_Uniform(self):
        for i in range(1):
            print("Uniform Test ",i)
            self.A, self.B, self.masses_A, self.masses_B = generate_points(N,DIM,"Uniform")
            distance_function = utils.dist
            self.hungarian = Hungarian(self.A, self.B, distance_function, P)
            self.hungarian.reset_matching()
            self.hungarian.compute_min_cost_matching()
            dist_matrix = [ [ pow(utils.dist(a,b),P) for a in self.A] for b in self.B ] 
            real_cost = pow(ot.emd2(self.masses_A, self.masses_B, dist_matrix), 1.0/P)
            cost = self.hungarian.min_cost
            ratio = cost / real_cost
            print()
            print("our matching cost: ", cost)
            print("emd2 matching cost: ", real_cost)
            print("ratio: ", ratio)

