import unittest
from constants import *
from distribution import generate_points
import utils
from hungarian import Hungarian
#from bottleneck import BottleneckMatching
import ot

class TestBottleneck(unittest.TestCase):

    def test_min_cost_matching(self):
        #self.p = P

        # GENERATE THE POINT SET AND DECOMPOSITION
        self.A, self.B, self.masses_A, self.masses_B = generate_points(400,10,"Uniform")
        distance_function = utils.l_1
        self.n = len(self.A)
        self.M = Hungarian(self.A, self.B, distance_function)
        self.diam = 2*pow(2, 10)

        #self.bottleneck.reset_matching()
        self.M.compute_min_cost_bottleneck_matching()
        #dist_matrix = [ [ pow(utils.l_infty(a,b),P) for a in self.A] for b in self.B ] 
        #self.M.emd_bottleneck()
        #real_cost = pow(ot.emd2(self.masses_A, self.masses_B, dist_matrix), 1.0)
        cost = self.M.min_cost
        #ratio = cost / real_cost
        print()
        print("our matching cost: ", cost)
        #print("emd2 matching cost: ", real_cost)
        #print("ratio: ", ratio)

#    '''
#    For testing
#    '''
#    def test_bottleneck(self):
        lb_cost, ub_cost = 0, self.diam
        bcount = 0
        while (ub_cost - lb_cost) > 0.01 and bcount < 1000:
            bcount += 1
            threshhold = (lb_cost + ub_cost)/2
            dist_matrix = [ [ 1 if utils.dist(a,b) <= threshhold else 2*self.n for a in self.A] for b in self.B ] 
            tcost = ot.emd2(self.masses_A, self.masses_B, dist_matrix)
            if tcost > self.n + 1:
                lb_cost = threshhold
            else:
                ub_cost = threshhold
        #print(threshhold)
        #print("cost", cost)
        print("emd matching cost: ", threshhold)
        ratio = cost / threshhold 
        print("ratio: ", ratio)
