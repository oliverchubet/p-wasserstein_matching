from constants import *
from distribution import generate_points
from utils import dist
from decomp import Decomposition
from itertools import chain
from tree import Tree
import copy
from math import ceil, log

class Hungarian:
    def __init__(self, A, B, dist, p=1, delta=0.01, base=1.01, dim=2):
        self.A = A
        self.B = B
        self.n = len(A)
        self.dist = dist
        self.p = p
        self.diam = 2*pow(2, p)
        self.min_cost = 0
        self.cost_using_clustering = 0
        self.decomp = Decomposition(A, B, dist, p, delta, k=2, base=base, dim=dim)
        self.visitedA = set()
        self.visitedB = set()
        self.freeB = set([b for b in self.B])
        self.unvisitedA = set([a for a in self.A])
        self.matched = 0
        self.phases = 0
        self.edges_seen = 0
        self.delta = delta
        self.base = base

    '''
    Computes min cost matching using p-th power distances
    '''
    def compute_min_cost_matching(self):
        if VERBOSE: print("\nCompute min-cost matching")
        count, len_path = 0, 0
        while self.matched < self.n and count < self.n: 
            count += 1
            print("Iteration", count)
            self.compute_cluster_distance()
            self.all_feasible()
            self.update_dual_weights(len_path)
            len_path = self.djikstra()
        if not self.is_matching():
            # Error message
            print("Algorithm finished without finding a perfect matching")
        self.compute_cluster_distance()

        self.all_feasible()
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
                b.dual_weight += shortest_path - b.len_path
            for a in self.visitedA:
                a.dual_weight += a.len_path - shortest_path
        self.unvisitedA = set([a for a in self.A])
        self.visitedA = set()
        self.visitedB = set([b for b in self.freeB])

    def all_feasible(self):
        for a in self.A:
            for b in self.B:
                assert(self.is_feasible(a,b))

    def is_feasible(self, a, b):
        proxy_dist = ceil(pow(self.distC[(a,b)],self.p)/self.delta)
        if a.dual_weight + b.dual_weight > proxy_dist + 1:
            print("Not feasible...")
            print("proxy dist= ", proxy_dist)
            print("y(a) = ", a.dual_weight)
            print("y(b) = ", b.dual_weight)
            return False

        if a.match is b:
            if a.dual_weight + b.dual_weight != proxy_dist:
                print("Not admissible...")
                print("proxy dist= ", proxy_dist)
                print("y(a) = ", a.dual_weight)
                print("y(b) = ", b.dual_weight)
                return False
        return True

    '''
    Shortest path tree used in Djikstra's algorithm
    '''
    def init_shortest_path_tree(self):
        self.shortest_path_tree = Tree()
        for b in self.freeB:
            self.shortest_path_tree.connect(b, "root")

    '''
    Compute the shortest path in the residual graph using slacks for edge costs
    Return an augmenting path and length
    '''
    def djikstra(self, threshhold=None):
        self.init_shortest_path_tree()
        self.decomp.init_min_slack_heap(self.unvisitedA, self.freeB, threshhold=threshhold)
        shortest_path = None
        while len(self.decomp.slack_heap) > 0:
            self.edges_seen += 1
            edge = self.decomp.find_min_slack_edge()
            a, b = edge.a, edge.b
            a.len_path = edge.slack + b.len_path
            self.shortest_path_tree.connect(a, b)
            self.visitedA.add(a)
            self.visitedB.add(b)
            self.unvisitedA.remove(a)
            if a.match is None:
                aug_path = self.shortest_path_tree.path(a)
                shortest_path = a.len_path
                #self.decomp.removeA(a, threshhold=threshhold)
                self.augment(aug_path)
                return shortest_path
            else:
                b = a.match
                b.len_path = a.len_path
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
        while a != "root":
            self.phases += 1
            b = next(it)
            b.dual_weight -= 1
            self.decomp.removeB(b)
            a.match = b
            a = next(it)
            self.shortest_path_tree.remove(a)
            self.shortest_path_tree.remove(b)
        if b in self.freeB:
            self.freeB.remove(b)
        self.matched += 1

    '''
    Used in testing
    '''
    def reset_matching(self):
        self.matched = 0
        for pt in self.A:
            pt.dual_weight = 0
            pt.len_path = 0
            pt.len_path = None
            pt.match = None
        for pt in self.B:
            pt.dual_weight = 0
            pt.len_path = 0
            pt.len_path = 0
            pt.match = None
            self.visitedB.add(pt)

    def compute_cluster_distance(self):
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


if __name__ == "__main__":
    A, B, masses_A, masses_B = generate_points(N,DIM,DISTRIBUTION)
    hungarian = Hungarian(A, B, dist, P)
