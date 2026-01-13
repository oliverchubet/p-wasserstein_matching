from constants import *
from utils import dist
from decomp import Decomposition
from itertools import chain
from tree import Tree
import copy
from math import ceil, log

class Wasserstein:
    def __init__(self, A, B, dist, p=1, delta=0.01, base=1.01):
        self.A = A
        self.B = B
        self.n = len(A)
        self.dist = dist
        self.p = p
        self.diam = 2*pow(2, DIM)
        self.min_cost = 0
        self.cost_using_clustering = 0
        self.decomp = Decomposition(A, B, dist, p, delta, k=2)
        self.freeB = set([b for b in self.B])
        self.matched = 0
        self.delta = delta

    def compute_pWasserstein(self):
        count = 0
        while self.matched < self.n and count < self.n:
            count += 1
            print("Iteration", count)
            self.dual_adjustment()
            self.augmentation()
            assert(self.all_feasible())
        assert(self.matched == self.n)
        self.compute_cluster_distance()
        cost = sum([pow(self.distC[(a,a.match)],self.p) for a in self.A])
        self.cost_using_clustering = pow(cost, 1.0/self.p)
        cost = sum([pow(dist(a,a.match),self.p) for a in self.A])
        self.min_cost = pow(cost, 1.0/self.p)

    def dual_adjustment(self):
        self.unvisitedA = set([a for a in self.A])
        self.visitedA = set()
        self.visitedB = set([a for a in self.freeB])
        for b in self.B: b.len_path = 0
        for a in self.A: a.len_path = None
        self.decomp.init_min_slack_heap(self.unvisitedA, self.freeB)
        if DEBUG: self.print_report()
        shortest_path = None
        count = 0
        edge = self.decomp.weighted_BCP()
        while edge is not None and shortest_path is None and count < self.n:
            count += 1
            a, b = edge.a, edge.b
            if DEBUG: print("edge", a.id, b.id, "slack", edge.slack)
            a.len_path = edge.slack
            self.visitedA.add(a)
            self.visitedB.add(b)
            self.decomp.removeA(a)
            if DEBUG: self.print_report()
            if a.match is not None:
                b = a.match
                self.visitedB.add(b)
                b.len_path = a.len_path
                self.decomp.updateB(b)
            else:
                shortest_path = a.len_path
            edge = self.decomp.weighted_BCP()
        if shortest_path:
            for a in self.visitedA:
                if DEBUG:
                    print("visitedA:")
                    print("y(", a.id,") =", a.dual_weight, "+", a.len_path, "-", shortest_path)
                a.dual_weight = a.dual_weight - shortest_path + a.len_path
            for b in self.visitedB:
                if DEBUG:
                    print("visitedB:")
                    print("y(",b.id,") =", b.dual_weight, "+", shortest_path, "-", b.len_path)
                b.dual_weight = b.dual_weight + shortest_path - b.len_path
        assert(self.all_feasible())

    def augmentation(self):
        if DEBUG: matched_at_least_one = False
        paths = []
        for b in self.freeB:
            path = [b]
            edge = self.decomp.weighted_NN(b)
            if DEBUG: print(edge)
            while edge is not None and edge.slack == 0 and edge.a.match is not None:
                path.extend([edge.a,edge.a.match])
                edge = self.decomp.weighted_NN(a.match)
            paths.append(path)
        for path in paths:
            if len(path) > 1 and path[-1].match is None:
                if DEBUG: matched_at_least_one = True
                it = iter(path)
                a = next(it)
                if DEBUG: print("augment:", end="")
                while a != "root":
                    b = next(it)
                    if DEBUG: print("-", a.id, "-", b.id, end=" ")
                    b.dual_weight -= 1
                    a.match = b
                    a = next(it)
                    self.shortest_path_tree.remove(a)
                    self.shortest_path_tree.remove(b)
                if DEBUG: print()
                self.freeB.remove(b)
                self.matched += 1
        if DEBUG: assert(matched_at_least_one)

    def all_feasible(self):
        self.compute_cluster_distance()
        #for a in self.A:
        #    print("a.id",a.id, "y(a)", a.dual_weight)
        #for b in self.B:
        #    print("b.id", b.id, "y(b)", b.dual_weight)
        #for a in self.A:
        #    if a.match is not None:
        #        print(a.id, "matched to", a.match.id)
        for a in self.A:
            for b in self.B:
                assert(self.is_feasible(a,b))
        return True

    def is_feasible(self, a, b):
        proxy_dist = ceil(pow(self.distC[(a,b)],self.p)/self.delta)
        check = True
        if a.dual_weight + b.dual_weight > proxy_dist + 1:
            print("Not feasible...")
            print("proxy dist= ", proxy_dist)
            print("y(",a.id, ") = ", a.dual_weight)
            print("y(", b.id, ") = ", b.dual_weight)
            print("y(a) + y(b) = ", a.dual_weight + b.dual_weight)
            index = self.distmap[(a,b)]
            print("distmap id", index)
            
            if a in self.unvisitedA: print("a unvisited")
            else: print("a visited")
            if b in self.visitedB: print("b visited")
            else: print("b unvisited")
            check = False
        if a.match is b:
            if a.dual_weight + b.dual_weight != proxy_dist:
                print("Not admissible...")
                print("proxy dist= ", proxy_dist)
                print("a.id", a.id, "y(a) = ", a.dual_weight)
                print("b.id", b.id, "y(b) = ", b.dual_weight)
                print("distmap id", self.distmap[(a,b)])
                check = False
        return check

    def compute_cluster_distance(self):
        # COMPUTE CLUSTER DISTANCE TO USE FOR MATCHING COST
        self.distC = dict()
        self.distmap = dict()
        for a in self.A:
            for b in self.B:
                self.distC[(a,b)] = self.diam 

        # BRUTE FORCE COMPUTE CLUSTER DISTANCE
        for center, cluster in self.decomp.clusters.items():
            for a in cluster.A:
                for b in cluster.B:
                    if self.distC[(a,b)] > cluster.distC(a,b):
                        self.distC[(a,b)] = cluster.distC(a,b)
                        self.distmap[(a,b)] = cluster.center.id
        #cost = sum([pow(distC[(a,a.match)],self.p) for a in self.A])

    def print_report(self):
        #for a in self.A:
            #print("y(",a.id,") =", a.dual_weight)
            #if a.match is not None:
            #    print("(",a.id,a.match.id,") in M")
        #for b in self.B:
            #print("y(",b.id,") =", b.dual_weight)

        for cluster in self.decomp.slack_heap._itemmap.keys():
            print("cluster", cluster.center.id, "edge", cluster.min_slack_edge.a.id, cluster.min_slack_edge.b.id, "min slack", cluster.min_slack_edge.slack)
