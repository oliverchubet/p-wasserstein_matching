import unittest
from distribution import generate_points
import utils
from random import randint
from wasserstein import Wasserstein
#from hungarian import Hungarian
from point import Point
from math import ceil

class TestMinimal(unittest.TestCase):

    def setUp(self):
        n = 4
        self.base = 1.01
        #self.delta = 1/(10*n)
        self.delta = 0.001
        self.p = 2
        self.A, self.B, masses_A, masses_B = generate_points(n,self.p,"Uniform")
        self.distance_function = utils.dist
        #for a in self.A:
        #    print(a.id,":", a.coords)
        #for b in self.B:
        #    print(b.id, ":", b.coords)
        self.wasserstein = Wasserstein(self.A, self.B, self.distance_function, self.p, delta=self.delta, base=self.base)
        self.wasserstein.compute_pWasserstein()
        #for a in self.A:
        #    for b in self.B:
        #        print(a.id, b.id, ceil(pow(self.wasserstein.distC[(a,b)],self.p)/self.delta))

    #def test_cluster_dist(self):
    #    self.wasserstein.compute_cluster_distance()
    #    for a in self.A:
    #        for b in self.B:
    #            dist = self.distance_function(a,b)
    #            distC = self.wasserstein.distC[(a,b)]
    #            assert(dist <= distC)
    #            assert(distC <= 4*self.base*dist)
                
    def test_wasserstein(self):
        self.wasserstein.compute_pWasserstein()

