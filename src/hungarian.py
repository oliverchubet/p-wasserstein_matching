from constants import *
from distribution import generate_points
from utils import dist
from decomp import Decomposition
from itertools import chain
from tree import Tree
import copy
from math import ceil, log

class Hungarian:
    def __init__(self, A, B, dist, p=1, delta=0.01):
        self.A = A
        self.B = B
        self.n = len(A)
        self.dist = dist
        self.p = p
        self.diam = 2*pow(2, DIM)
        self.min_cost = 0
        self.cost_using_clustering = 0
        self.decomp = Decomposition(A, B, dist)
        self.visitedA = set()
        self.visitedB = set()
        self.freeB = set([b for b in self.B])
        self.matched = 0
        self.phases = 0
        self.edges_seen = 0
        self.delta = delta
        self.danger = 0

    '''
    Computes min cost matching using p-th power distances
    '''
    def compute_min_cost_matching(self):
        if VERBOSE: print("\nCompute min-cost matching")
        count = 0
        while self.matched < self.n and count < self.n: 
            count += 1
            if VERBOSE: print("x", end='', flush=True)
            len_path = self.djikstra()
            self.update_dual_weights(len_path)
            if self.danger > 5:
                print("Holy fuck")
        if not self.is_matching():
            # Error message
            print("Algorithm finished without finding a perfect matching")
        # COMPUTE CLUSTER DISTANCE TO USE FOR MATCHING COST
        self.distC = dict()
        for a in self.A:
            for b in self.B:
                self.distC[(a,b)] = self.diam 

        # BRUTE FORCE COMPUTE CLUSTER DISTANCE
        for center, cluster in self.decomp.clusters.items():
            for a in cluster.A:
                for b in cluster.B:
                    self.distC[(a,b)] = min(cluster.distC(a,b), self.distC[(a,b)])
        #cost = sum([pow(distC[(a,a.match)],self.p) for a in self.A])

        self.check_all_feasible()
        cost = sum([pow(self.distC[(a,a.match)],self.p) for a in self.A])
        self.cost_using_clustering = pow(cost, 1.0/self.p)
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
                ub_cost = max([self.dist(a,a.match) for a in self.A])
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
        if shortest_path is not None:
            for b in self.visitedB:
                if b.len_path < shortest_path:
                        b.dual_weight += shortest_path - b.len_path
                        b.weight = b.dual_weight
            for a in self.visitedA:
                if a.len_path < shortest_path:
                    a.dual_weight += a.len_path #- shortest_path
                    a.weight = a.dual_weight
        else:
            self.danger += 1
        self.visitedB = set([b for b in self.freeB])
        self.visitedA = set()

    def check_all_feasible(self):
        for a in self.A:
            for b in self.B:
                self.check_feasible(a,b)

    def check_feasible(self, a, b):
        level = ceil(log(self.dist(a,b),BASE))

        slack = pow(2*pow(BASE, level),self.p) - a.weight - b.weight + self.delta
        if -a.dual_weight + b.dual_weight > slack:
            print("Not feasible...")
            print("slack = ", slack)
            print("y(a) = ", a.dual_weight)
            print("y(b) = ", b.dual_weight)

        if a.match is b:
            if -a.dual_weight + b.dual_weight > slack:
                print("Not admissible...")
                print("slack = ", slack)
                print("y(a) = ", a.dual_weight)
                print("y(b) = ", b.dual_weight)

    '''
    Shortest path tree used in Djikstra's algorithm
    '''
    def init_shortest_path_tree(self):
        self.shortest_path_tree = Tree()
        for b in self.visitedB:
            self.shortest_path_tree.connect(b, "root")

    '''
    Compute the shortest path in the residual graph using slacks for edge costs
    Return an augmenting path and length
    '''
    def djikstra(self, threshhold=None):
        self.init_shortest_path_tree()
        self.decomp.init_min_slack_heap(self.visitedA, self.visitedB, threshhold=threshhold)
        len_path = self.diam 
        while len(self.decomp.slack_heap) > 0:
            #if VERBOSE: print(".", end='', flush=True)
            self.edges_seen += 1
            edge = self.decomp.find_min_slack_edge()
            if edge.slack is None:
                return len_path
            a, b = edge.a, edge.b
            a.len_path = edge.slack + b.len_path
            self.shortest_path_tree.connect(a, b)
            self.visitedA.add(a)
            if a.free:
                #if VERBOSE: print("x",self.dist(a,b), end='', flush=True)
                aug_path = self.shortest_path_tree.path(a)
                len_path = a.len_path
                self.decomp.removeA(a, threshhold=threshhold)
                self.augment(aug_path)
                #return a.len_path
            else:
                b = a.match
                b.len_path = a.len_path
                b.weight = b.dual_weight - b.len_path
                self.shortest_path_tree.connect(b, a)
                self.visitedB.add(b)
                self.decomp.removeA(a, threshhold=threshhold)
                self.decomp.updateB(b, threshhold=threshhold)
        #return len_path

    '''
    Updates the matching given the path returned by Djikstra's
    '''
    def augment(self, path):
        it = iter(path)
        a = next(it)
        a.free = False
        while a != "root":
            self.phases += 1
            b = next(it)
            b.dual_weight -= self.delta
            b.weight -= self.delta
            self.decomp.removeB(b)
            a.match = b
            a = next(it)
            self.shortest_path_tree.remove(a)
            self.shortest_path_tree.remove(b)
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
