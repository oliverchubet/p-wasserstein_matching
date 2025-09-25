from constants import *
from distribution import generate_points
from utils import dist
from decomp import Decomposition
from itertools import chain
from tree import Tree
import copy

class Hungarian:
    def __init__(self, A, B, dist, p=1):
        self.A = A
        self.B = B
        self.n = len(A)
        self.dist = dist
        self.p = p
        self.diam = 2*pow(2, DIM)
        self.min_cost = 0
        self.decomp = Decomposition(A, B, dist)
        self.visitedA = set()
        self.visitedB = set([b for b in self.B])
        self.freeB = set([b for b in self.B])
        self.matched = 0
        self.ops = 0

    '''
    Computes min cost matching using p-th power distances
    '''
    def compute_min_cost_matching(self):
        if VERBOSE: print("\nCompute min-cost matching")
        count = 0
        while self.matched < self.n and count < self.n: 
            count += 1
            if VERBOSE: print("\n", count, end='', flush=True)
            len_path = self.djikstra()
            self.update_dual_weights(len_path)
        if not self.is_matching():
            # Error message
            print("Algorithm finished without finding a perfect matching")
        cost = sum([pow(dist(a,a.match),self.p) for a in self.A])
        self.min_cost = pow(cost, 1.0/self.p)

    def compute_min_cost_bottleneck_matching(self):
        if VERBOSE: print("\nCompute min-cost bottleneck matching")
        lb_cost, ub_cost = 0, self.diam
        count = 0
        bcount = 0
        while (ub_cost - lb_cost) > 0.01 and bcount < 1000:
            self.reset_matching()
            if VERBOSE: print("\nLB cost =", lb_cost, "UB cost =", ub_cost)
            cost = (lb_cost + ub_cost)/2
            count = 0
            bcount += 1
            while self.matched < self.n and count < self.n: 
                if self.matched < count: break
                count += 1
                if VERBOSE: print("\n", count, end='', flush=True)
                len_path = self.djikstra(threshhold=cost)
                self.update_dual_weights(len_path)
            if not self.is_matching():
                lb_cost = cost
            else:
                ub_cost = max([dist(a,a.match) for a in self.A])
            self.min_cost = ub_cost

    def is_matching(self):
        for a in self.A:
            if a.match is None:
                return False
        return True

    '''
    Uses the shortest path length returned by Djikstra's to update the dual weights
    '''
    def update_dual_weights(self, shortest_path):
        for b in self.visitedB:
            if b.len_path < shortest_path:
                    b.dual_weight += shortest_path - b.len_path
                    b.weight = b.dual_weight - b.len_path
        for a in self.visitedA:
            if a.len_path < shortest_path:
                a.dual_weight += a.len_path - shortest_path
        for b in self.visitedB: b.in_tree = False
        self.visitedB = set([b for b in self.freeB])
        self.visitedA = set()

    '''
    Shortest path tree used in Djikstra's algorithm
    '''
    def init_shortest_path_tree(self):
        self.shortest_path_tree = Tree()
        for b in self.visitedB:
            self.shortest_path_tree.connect(b, "root")
            b.in_tree = True

    '''
    Compute the shortest path in the residual graph using slacks for edge costs
    Return an augmenting path and length
    '''
    def djikstra(self, threshhold=None):
        self.init_shortest_path_tree()
        self.decomp.init_min_slack_heap(self.visitedA, self.visitedB, threshhold=threshhold)
        len_path = self.diam 
        while len(self.decomp.slack_heap) > 0:
            if VERBOSE: print(".", end='', flush=True)
            edge = self.decomp.find_min_slack_edge()
            a, b = edge.a, edge.b
            a.len_path = edge.slack + b.len_path
            self.shortest_path_tree.connect(a, b)
            self.visitedA.add(a)
            if a.free:
                if VERBOSE: print("x",self.dist(a,b), end='', flush=True)
                aug_path = self.shortest_path_tree.path(a)
                len_path = a.len_path
                self.decomp.removeA(a, threshhold=threshhold)
                self.augment(aug_path)
                return a.len_path
            else:
                b = a.match
                b.len_path = a.len_path
                b.weight = b.dual_weight - b.len_path
                self.shortest_path_tree.connect(b, a)
                self.visitedB.add(b)
                self.decomp.removeA(a, threshhold=threshhold)
                self.decomp.updateB(b, threshhold=threshhold)
        return len_path

    '''
    Updates the matching given the path returned by Djikstra's
    '''
    def augment(self, path):
        it = iter(path)
        a = next(it)
        a.free = False
        while a != "root":
            self.ops += 1
            b = next(it)
            a.match = b
            a = next(it)
        if b in self.freeB:
            self.freeB.remove(b)
            b.free = False
        self.matched += 1

    '''
    Used in testing
    '''
    def reset_matching(self):
        self.matched = 0
        self.visitedA = set()
        self.visitedB = set()
        for pt in self.A:
            pt.dual_weight = 0
            pt.weight = 0
            pt.len_path = None
            pt.match = None
            pt.visited = False
            pt.in_tree = False
            pt.free = True
        for pt in self.B:
            pt.dual_weight = 0
            pt.weight = 0
            pt.len_path = 0
            pt.match = None
            self.visitedB.add(pt)



if __name__ == "__main__":
    A, B, masses_A, masses_B = generate_points(N,DIM,DISTRIBUTION)
    hungarian = Hungarian(A, B, dist, P)
