import unittest
import matplotlib.pyplot as plt
from decomp import Decomposition
from collections import defaultdict
from itertools import chain
from math import log, ceil
import ot
from constants import *
from distribution import generate_points
import math
import utils

# TESTS RELATED TO THE CLUSTER DECOMPOSITION
class TestClusterParticipation(unittest.TestCase):

    def test_cluster_participation(self):

        dim = 2
        diam = 2*pow(2, 10)
        dist = utils.dist

        x_avg_k2 = list()
        y_avg_k2 = list()

        #for n in [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000]:
        for n in [20]:
            #for dim in range(1,11):
            for p in [2]:
                for distr in ['Uniform']:
                    print("n =", n, "dim = ", dim, "p = ", p, distr)
                    A, B, masses_A, masses_B = generate_points(n,dim,distr)
                    decomp = Decomposition(A, B, dist, p=p, k=2)
                    data = []
                    # INITIALIZE DICTIONARY FOR CLUSTER DISTANCE
                    distC = dict()
                    for a in A:
                        for b in B:
                            distC[(a,b)] = diam 

                    # BRUTE FORCE COMPUTE CLUSTER DISTANCE
                    for center, cluster in decomp.clusters.items():
                        for a in cluster.A:
                            for b in cluster.B:
                                distC[(a,b)] = min(cluster.distC(a,b), distC[(a,b)])
                                
                    # CHECK APPROXIMATION FACTOR OF CLUSTER DISTANCE
                    max_approx = 1
                    count = 0
                    total_approx = 0
                    for a in A:
                        for b in B:
                            count += 1
                            approx = distC[(a,b)]/utils.dist(a,b)

                            max_approx = max(max_approx, approx)
                            total_approx += approx
                    avg_approx = total_approx / count

                    entry = "(" + ", ".join([str(2*n), str(dim), str(p), str(distr), str(max_approx), str(avg_approx)]) + "),\n"
                    print("Result:", entry)
                    with open("distortion_d2_l2.py", "a") as myfile:
                        myfile.write(entry)
