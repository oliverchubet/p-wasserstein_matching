import unittest
from distribution import generate_points
from heaps import MaxLevelHeap
from collections import defaultdict

class TestHeap(unittest.TestCase):

    def test_heap(self):
        self.p = 2
        self.delta = 0.01
        self.k = 2
        self.n = 5
        self.dim = 2
        self.base = 1.01
        self.distribution = "Normal"
        self.A, self.B, self.masses_A, self.masses_B = generate_points(self.n,self.dim,self.distribution)
        self.heapB = MaxLevelHeap()
        for b in self.B:
            b.len_path = 0
        buckets = defaultdict(set)
        buckets[None] = {self.B[0]}
        buckets[-5] = {self.B[1], self.B[2], self.B[3]}
        buckets[-10] = {self.B[4]}
        for level, bucket in buckets.items():
            for b in bucket:
                #if level not in self.heapB.heaps:
                    #print(f"{level} not in heaps")
                self.heapB.insert(b, b.dual_weight - b.len_path,level)
        #print(self.heapB)
        #print(self.heapB.findmax().id)
        #print("max_val", self.heapB.max_val)
        #print("max_item", self.heapB.max_item.id)

        b = self.B[0]
        self.heapB.changepriority(b, -10)
        b = self.B[1]
        self.heapB.changepriority(b, -10)
        b = self.B[3]
        self.heapB.changepriority(b, -10)
        b = self.B[4]
        self.heapB.changepriority(b, -10)
        #print(self.heapB.findmax().id)
        #print(self.heapB)
        #print("max_val", self.heapB.max_val)
        #print("max_item", self.heapB.max_item.id)

