from constants import *
from itertools import chain
from collections import defaultdict
from math import sqrt, log, ceil
import random
import math
#from ds2.priorityqueue import PriorityQueue as MinHeap
from heaps import MaxLevelHeap, MinHeap
from point import Edge
import utils

class Decomposition:
    def __init__(self, A, B, dist, p, delta, k, base=BASE, dim=DIM):
        if VERBOSE: print("Compute decomposition")
        self.A = A
        self.B = B
        self.n = len(A) + len(B)
        self.dist = dist
        self.p = p
        self.delta = delta
        self.k = k
        self.diam = pow(2, dim)
        self.base = base
        self.layers = defaultdict(set)
        self.clusters = dict()
        self.lookup = defaultdict(set)
        self.slack_heap = None
        self.init_nearest()

        self.compute_samples()
        for i in range(self.k-1, 0, -1):
            if VERBOSE: print("Computing layer", i, " clusters", end=" ")
            self.compute_clusters(i)
            self.compute_nearest(i)
        if VERBOSE: print("\nComputing layer 0 clusters", end=" ")
        self.compute_clusters(0)
        if VERBOSE: print()

    '''
    Each point in A and B is assigned to a layer
    '''
    def compute_samples(self):
        prob = pow(self.n, -1.0/self.k)
        for pt in chain(self.A, self.B):
            for i in range(self.k):
                if i == self.k-1 or random.random() > prob:
                    self.layers[i].add(pt)
                    break

    '''
    Initializes self.nearest to store nearest neighbor distances to the current sample
    '''
    def init_nearest(self):
        self.nearest = dict()
        for x in chain(self.A,self.B):
            self.nearest[x] = self.diam

    '''
    Brute force computation of nearest neighbor distances to sample
    '''
    def compute_nearest(self, i):
        for x in chain(self.A,self.B):
            for s in self.layers[i]:
                self.nearest[x] = min(self.nearest[x], self.dist(x,s))

    def compute_clusters(self, i):
        for c in self.layers[i]:
            self.compute_cluster(c)

    '''
    Adds points to the cluster if the center is closer than the nearest neighbor distance to the sample.
    Points are bucketed by levels determined by distance to the center.
    '''
    def compute_cluster(self, center):
        if VERBOSE: print('.', end='', flush=True)
        bucketsA, bucketsB = defaultdict(set), defaultdict(set)
        pointsA, pointsB = dict(), dict()
        for a in self.A:
            if self.dist(a,center) < self.nearest[a]:
                self.lookup[a].add(center)
                level = None
                if self.dist(a,center) > 0:
                    level = ceil(log(self.dist(a,center),self.base))
                bucketsA[level].add(a)
                pointsA[a] = level
        for b in self.B:
            if self.dist(center,b) < self.nearest[b]:
                self.lookup[b].add(center)
                level = None
                if self.dist(center,b) > 0:
                    level = ceil(log(self.dist(center,b), self.base))
                bucketsB[level].add(b)
                pointsB[b] = level
        if len(bucketsA) > 0 and len(bucketsB) > 0:
            self.clusters[center] = Cluster(pointsA, pointsB,center, bucketsA, bucketsB, self.p, self.base, self.delta, self.diam)

    '''
    Take the min slack edge from each cluster and put it into a min heap
    '''
    def init_min_slack_heap(self, A, B): #threshhold=None):
        slack_heap = MinHeap()
        for center, cluster in self.clusters.items():
            cluster.init_max_heaps(A, B)
            #cluster.print_heap_contents()
            edge = cluster.weighted_BCP()
            if edge is not None:
                slack_heap.insert(cluster, edge.slack)
        self.slack_heap = slack_heap


    '''
    Returns the min slack edge but does not remove it from the heap
    '''
    def weighted_BCP(self):
        if len(self.slack_heap) > 0:
            cluster = self.slack_heap.findmin()
            if DEBUG:
                print("min slack edge center id", cluster.center.id)
            return cluster.min_slack_edge

    def weighted_NN(self, x):
        min_slack_edge = None
        min_slack = self.diam
        for center in self.lookup[x]:
            if center not in self.clusters:
                continue
            edge = self.clusters[center].weighted_NN(x)
            if edge is None:
                continue
            if min_slack_edge is None or edge.slack < min_slack:
                min_slack_edge = edge
                min_slack = edge.slack
        return min_slack_edge

    '''
    Each loop of Djikstra's, some points will have their path length updated.
    Tell the clusters affected to update their heaps and min slack edge
    Updating the min slack heap
    '''
    def removeA(self, a): #threshhold=None):
        for center in self.lookup[a]:
            if center not in self.clusters:
                continue
            cluster = self.clusters[center]
            if cluster not in self.slack_heap:
                continue
            cluster.removeA(a)
            edge = cluster.weighted_BCP()
            if edge is None:
                self.slack_heap.remove(cluster)
            else:
                self.slack_heap.changepriority(cluster, edge.slack)

    def removeB(self, b): #threshhold=None):
        for center in self.lookup[b]:
            if center not in self.clusters:
                continue
            cluster = self.clusters[center]
            if cluster not in self.slack_heap:
                continue
            cluster.removeB(b)
            edge = cluster.weighted_BCP() #threshhold=threshhold)
            if edge is None:
                self.slack_heap.remove(cluster)
            else:
                self.slack_heap.changepriority(cluster, edge.slack)

    def updateB(self, b): #threshhold=None):
        for center in self.lookup[b]:
            if center not in self.clusters:
                continue 
            cluster = self.clusters[center]
            if cluster not in self.slack_heap:
                continue
            cluster.updateB(b)
            edge = cluster.weighted_BCP() #threshhold=threshhold)
            if edge is None:
                self.slack_heap.remove(cluster)
            else:
                self.slack_heap.changepriority(cluster, edge.slack)

    def print_clusters(self):
        print("--------------------------------")
        for cluster in self.clusters:
            self.clusters[cluster].print_cluster()
        print("--------------------------------")


class Cluster:
    def __init__(self, A, B, center, bucketsA, bucketsB, p, base, delta, diam):
        self.A = A
        self.B = B
        self.center = center
        self.bucketsA = bucketsA
        self.bucketsB = bucketsB
        self.p = p
        self.base = base
        self.delta = delta
        self.diam = diam

     '''
     Once the max heaps and max weights are initialized, the min slack edge is computed
     '''
    def weighted_BCP(self):
        max_a = None
        max_b = None
        min_slack = pow(2, self.p)
#+        self.print_heap_contents()
#+        max_a, max_b = None, None
#+        min_slack = None
         min_slack_edge = None
        a,b = None, None
        if threshhold is not None: tlevel = ceil(log(threshhold,self.base))
         for level in self.levels:
            if threshhold is None or level <= tlevel:
                if level in self.heapsA and len(self.heapsA[level]) > 0:
                    a = self.heapsA[level].findmax()
                if level in self.heapsB and len(self.heapsB[level]) > 0:
                    b = self.heapsB[level].findmax()
                if max_a is None or a.weight > max_a.weight:
                    max_a = a
                if max_b is None or b.weight > max_b.weight :
                    max_b = b
                if a is not None and b is not None:
                    slack = pow(2*pow(self.base, level),self.p) - max_a.weight - max_b.weight
                    if min_slack_edge is None or min_slack_edge.slack > slack:
                        min_slack_edge = Edge(max_a, max_b, slack, self.center, level)
#+            print(level)
#+            max_a = self.heapA.findmax(level)
#+            max_b = self.heapB.findmax(level)
#+            if max_a is not None and max_b is not None:
#+                slack = self.proxyDistC(max_a, max_b) - max_a.dual_weight - (max_b.dual_weight - max_b.len_path) + 1
#+                if min_slack is None or min_slack > slack:
#+                    min_slack = slack
#+                    min_slack_edge = Edge(max_a, max_b, slack, self.center, level)
#+                    print(max_a.id, max_b.id)
         self.min_slack_edge = min_slack_edge
         return min_slack_edge
#    def weighted_BCP(self):
#        self.print_heap_contents()
#        max_a, max_b = None, None
#        min_slack = None
#        min_slack_edge = None
#        for level in self.levels:
#            print(level)
#            max_a = self.heapA.findmax(level)
#            max_b = self.heapB.findmax(level)
#            if max_a is not None and max_b is not None:
#                slack = self.proxyDistC(max_a, max_b) - max_a.dual_weight - (max_b.dual_weight - max_b.len_path) + 1
#                if min_slack is None or min_slack > slack:
#                    min_slack = slack
#                    min_slack_edge = Edge(max_a, max_b, slack, self.center, level)
#                    print(max_a.id, max_b.id)
#        self.min_slack_edge = min_slack_edge
#        return min_slack_edge

    def weighted_NN(self,b):
        max_a = None
        min_slack = self.diam
        min_slack_edge = None
        for level in self.levels:
            max_a = self.heapA.findmax(level)
            if max_a is not None:
                slack = self.proxyDistC(max_a, b) - max_a.dual_weight - (b.dual_weight - b.len_path)
                if min_slack > slack:
                    min_slack = slack
                    min_slack_edge = Edge(max_a, b, slack, self.center, level)
        return min_slack_edge

    def removeA(self, a):
        self.heapA.remove(a)

    def removeB(self, b):
        self.heapB.remove(b)

    def updateB(self, b):
        level = self.B[b]
        self.heapB.insert(b, b.dual_weight - b.len_path, level)
                
    '''
    Insert points into heaps by level keyed by their weight
    Initializes self.heapsA and self.heapsB (each level maps to a max heap keyed by weights)
    '''
    def init_max_heaps(self, A, B):
        self.heapA, self.heapB = MaxLevelHeap(), MaxLevelHeap()
        for level, bucket in self.bucketsA.items():
            for a in bucket:
                if a in A:
                    self.heapA.insert(a, a.dual_weight, level)
        for level, bucket in self.bucketsB.items():
            for b in bucket:
                if b in B:
                    self.heapB.insert(b, b.dual_weight - b.len_path,level)
        self.compute_levels()

    '''
    Computes a sorted list of the non-empty levels
    '''
    def compute_levels(self):
        self.levelsA = self.heapA.levels
        self.levelsB = self.heapB.levels
        levels = self.levelsA.union(self.levelsB)
        if None in levels: levels.remove(None)
        levels = list(levels)
        levels.sort()
        self.levels = [None] + levels

    '''
    Cluster distance within this cluster
    '''
    def distC(self, ptA, ptB):
        if self.A[ptA] is None : level = self.B[ptB]
        elif self.B[ptB] is None : level = self.A[ptA]
        else: level = max(self.A[ptA], self.B[ptB])
        if level is None: return 0
        return 2*pow(self.base, level)

    def proxyDistC(self, ptA, ptB):
        return ceil(pow(self.distC(ptA,ptB),self.p)/self.delta)
    def print_cluster(self):
        print("Cluster", self.center.id, end=": ")
        print(" ".join([str(a.id) for a in self.A]), end=" + ")
        print(" ".join([str(b.id) for b in self.B]))

    def print_heap_contents(self):
        print("Cluster", self.center.id, ":")
        print("Heap A:", self.heapA)
        print("Heap B:", self.heapB)
