from constants import *
from itertools import chain
from collections import defaultdict
from math import sqrt, log, ceil
import random
import math
#from ds2.priorityqueue import PriorityQueue as MinHeap
from heaps import MaxHeap, MinHeap
from point import Edge
import utils

class Decomposition:
    def __init__(self, A, B, dist, p=P, k=K, base=BASE):
        if VERBOSE: print("Compute decomposition")
        self.A = A
        self.B = B
        self.n = len(A) + len(B)
        self.dist = dist
        self.p = p
        self.k = k
        self.diam = pow(2, DIM)
        self.base = base
        self.layers = defaultdict(set)
        self.clusters = dict()
        self.lookup = defaultdict(set)
        self.init_nearest()

        self.compute_samples()
        for i in range(self.k-1, 0, -1):
            if VERBOSE: print("Computing layer", i, " clusters")
            self.compute_clusters(i)
            self.compute_nearest(i)
        if VERBOSE: print("Computing layer 0 clusters")
        self.compute_clusters(0)

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
            self.clusters[center] = Cluster(pointsA, pointsB,center, bucketsA, bucketsB, self.p, self.base)

    '''
    Take the min slack edge from each cluster and put it into a min heap
    '''
    def init_min_slack_heap(self, visitedA, visitedB, threshhold=None):
        slack_heap = MinHeap()
        for center, cluster in self.clusters.items():
            cluster.init_max_heaps(visitedA, visitedB)
            edge = cluster.compute_min_slack_edge(threshhold=threshhold)
            if edge is not None:
                slack_heap.insert(cluster, edge.slack)
        self.slack_heap = slack_heap


    '''
    Returns the min slack edge but does not remove it from the heap
    '''
    def find_min_slack_edge(self):
        if len(self.slack_heap) > 0:
            cluster = self.slack_heap.findmin()
            return cluster.min_slack_edge

    '''
    Each loop of Djikstra's, some points will have their path length updated.
    Tell the clusters affected to update their heaps and min slack edge
    Updating the min slack heap
    '''
    def removeA(self, a, threshhold=None):
        for center in self.lookup[a]:
            if center not in self.clusters:
                continue
            cluster = self.clusters[center]
            if cluster not in self.slack_heap:
                continue
            cluster.removeA(a)
            edge = cluster.compute_min_slack_edge(threshhold=threshhold)
            if edge is None:
                self.slack_heap.remove(cluster)
            else:
                self.slack_heap.changepriority(cluster, edge.slack)

    def removeB(self, b, threshhold=None):
        for center in self.lookup[b]:
            if center not in self.clusters:
                continue
            cluster = self.clusters[center]
            if cluster not in self.slack_heap:
                continue
            cluster.removeB(b)
            edge = cluster.compute_min_slack_edge(threshhold=threshhold)
            if edge is None:
                self.slack_heap.remove(cluster)
            else:
                self.slack_heap.changepriority(cluster, edge.slack)

    def updateB(self, b, threshhold=None):
        for center in self.lookup[b]:
            if center not in self.clusters:
                continue 
            cluster = self.clusters[center]
            if cluster not in self.slack_heap:
                continue
            cluster.updateB(b)
            edge = cluster.compute_min_slack_edge(threshhold=threshhold)
            if edge is None:
                self.slack_heap.remove(cluster)
            else:
                self.slack_heap.changepriority(cluster, edge.slack)


class Cluster:
    def __init__(self, A, B, center, bucketsA, bucketsB, p, base):
        self.A = A
        self.B = B
        self.center = center
        self.bucketsA = bucketsA
        self.bucketsB = bucketsB
        self.p = p
        self.base = base

    '''
    Once the max heaps and max weights are initialized, the min slack edge is computed
    '''
    def compute_min_slack_edge(self, threshhold=None):
        max_a = None
        max_b = None
        min_slack = pow(2, self.p)
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
        self.min_slack_edge = min_slack_edge
        return min_slack_edge

    def removeA(self, a):
        level = self.A[a]
        if a in self.heapsA[level]:
            self.heapsA[level].remove(a)

    def removeB(self, b):
        level = self.B[b]
        if b in self.heapsB[level]:
            self.heapsB[level].remove(b)

    def updateB(self, b):
        level = self.B[b]
        self.heapsB[level].insert(b, b.weight)
                
    '''
    Insert points into heaps by level keyed by their weight
    Initializes self.heapsA and self.heapsB (each level maps to a max heap keyed by weights)
    '''
    def init_max_heaps(self, visitedA, visitedB):
        self.heapsA, self.heapsB = dict(), dict()
        for level, bucket in self.bucketsA.items():
            heap = MaxHeap()
            for a in bucket:
                if not a in visitedA:
                    heap.insert(a, a.dual_weight)
            self.heapsA[level] = heap
        for level, bucket in self.bucketsB.items():
            heap = MaxHeap()
            for b in bucket:
                if b in visitedB:
                    heap.insert(b, b.dual_weight - b.len_path)
            self.heapsB[level] = heap
        self.compute_levels()

    '''
    Computes a sorted list of the non-empty levels
    '''
    def compute_levels(self):
        self.levelsA = set(self.heapsA.keys())
        self.levelsB = set(self.heapsB.keys())
        levels = self.levelsA.union(self.levelsB)
        if None in levels: levels.remove(None)
        levels = list(levels)
        levels.sort()
        self.levels = levels

    '''
    Used in testing
    '''
    def distC(self, ptA, ptB):
        if self.A[ptA] is None : level = self.B[ptB]
        elif self.B[ptB] is None : level = self.A[ptA]
        else: level = max(self.A[ptA], self.B[ptB])
        if level is None: return 0
        return 2*pow(self.base, level)

