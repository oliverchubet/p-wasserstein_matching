import unittest
from constants import *
#from distribution import generate_points
import utils
from hungarian import Hungarian
import ot
from tensorflow.keras.datasets import mnist
import numpy as np
from point import Point

class TestRealData(unittest.TestCase):

    def test_real_data_mnist(self):
        sample_sizes = [250, 750, 1500, 2000, 4000]
        dim = 728
        digits = [0,1,2,3,4,5,6,7,8,9]
        (x_train, y_train), (x_test, y_test) = mnist.load_data()
        for p in [2,5,10,20]:
            for n in sample_sizes:
                self.A, self.B = [], []
                for digit in digits:
                    train_mask = (y_train == digit)
                    x_train_digit = x_train[train_mask]
                    self.A.extend([Point(a,None) for a in x_train_digit[0:n//10]])
                    self.B.extend([Point(b,None) for b in x_train_digit[n//10:2*n//10]])
                self.masses_A = [1 for _ in range(n)]
                self.masses_B = [1 for _ in range(n)]
                distance_function = lambda a, b: np.linalg.norm(a.coords-b.coords)/256.0
                self.hungarian = Hungarian(self.A, self.B, distance_function, p)
                self.hungarian.reset_matching()
                #self.hungarian.p = p
                self.hungarian.compute_min_cost_matching()
                dist_matrix = [ [ pow(distance_function(a,b),p) for a in self.A] for b in self.B ] 
                real_cost = pow(ot.emd2(self.masses_A, self.masses_B, dist_matrix), 1.0/p)
                cost = self.hungarian.min_cost
                ratio = cost / real_cost
                print()
                print("our matching cost: ", cost)
                print("emd2 matching cost: ", real_cost)
                print("ratio: ", ratio)
                filename = "mnist_results"
                with open(filename, "a") as f:
                    f.write("\t".join([str(n), str(p), str(dim), str(real_cost), str(cost), str(ratio),"\n"]))


