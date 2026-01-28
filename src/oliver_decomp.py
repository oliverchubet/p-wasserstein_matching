from constants import *
from itertools import chain
from collections import defaultdict
from math import sqrt, log, ceil, inf
import random
import math
#from ds2.priorityqueue import PriorityQueue as MinHeap
from ai_heap import MaxLevelHeap 
from ai_heap import PriorityQueue as MinHeap
from point import Edge
import utils
from skipheap import SkipHeap

TRACE = False

class Decomposition:
    def __init__(self, A, B, dist, p, delta, k, base, dim):
        self.A = A
        self.B = B
        self.n = len(A) + len(B)
        self.dist = dist
        self.p = p
        self.delta = delta
        self.k = k
        self.diam = pow(2, dim)
        self.base = base
        self.layersA = defaultdict(set)
        self.layersB = defaultdict(set)
        self.levels = {-inf}
        self.clusters = dict()
        self.lookup = {x : defaultdict(set) for x in chain(A, B)}
        self.slack_heap = None
        self.init_nearest()
        self.inactivated_A = set()
        self.activated_B = set([b for b in self.B])

        self.compute_samples()
        for i in range(self.k-1, 0, -1):
            self.compute_clusters(i)
            self.compute_nearest(i)

        self.compute_clusters(0)

        self.cache_level_dist = {level: pow(base, level) for level in self.levels}
        self.cache_proxy_dist = {level: ceil(pow(self.cache_level_dist[level],p)/self.delta) for level in self.levels}
        self.cache_proxy_diam = {level: ceil(pow(2*self.cache_level_dist[level],p)/self.delta) for level in self.levels}
        self._max_slack = ceil(pow(self.diam, self.p)/self.delta)

        self.init_min_slack_heap()

    '''
    Each point in A and B is assigned to a layer
    '''
    def compute_samples(self):
        prob = pow(self.n, -1.0/self.k)
        #prob = log(self.n)/self.n

        for pt in chain(self.A):
            for i in range(self.k):
                self.layersA[i].add(pt)
                if i == self.k-1 or random.random() > prob:
                    break
        for pt in chain(self.B):
            for i in range(self.k):
                self.layersB[i].add(pt)
                if i == self.k-1 or random.random() > prob:
                    break


    '''
    Initializes self.nearest to store nearest neighbor distances to the current sample
    '''
    def init_nearest(self):
        self.nearestA, self.nearestB = dict(), dict()
        for x in chain(self.A,self.B):
            self.nearestA[x], self.nearestB[x] = self.diam, self.diam

    '''
    Brute force computation of nearest neighbor distances to sample
    '''
    def compute_nearest(self, i):
        for x in chain(self.A,self.B):
            for s in self.layersA[i]:
                self.nearestA[x] = min(self.nearestA[x], self.dist(x,s))
            for s in self.layersB[i]:
                self.nearestB[x] = min(self.nearestB[x], self.dist(x,s))

    def compute_clusters(self, i):
        for c in self.layersA[i]:
            self.compute_clusterA(c,i)
        for c in self.layersB[i]:
            self.compute_clusterB(c,i)

    '''
    Adds points to the cluster if the center is closer than the nearest neighbor distance to the sample.
    Points are bucketed by levels determined by distance to the center.
    '''
    def compute_clusterA(self, center, i):
        #if VERBOSE: print('.', end='', flush=True)
        bucketsA, bucketsB = defaultdict(set), defaultdict(set)
        pointsA, pointsB = dict(), dict()
        #for a in self.layersA[i]:
        for a in self.A:
            #if self.dist(a,center) < min(self.nearestA[a], self.nearestB[a]):
            #if (i == 1 and a in self.layersA[i]) or (i != 1 and self.dist(a,center) < self.nearestA[a]):
            if (i != 1 and self.dist(a,center) < self.nearestA[a]):
                level = -inf
                if self.dist(a,center) > 0:
                    level = ceil(log(self.dist(a,center),self.base))
                    self.levels.add(level)
                bucketsA[level].add(a)
                self.lookup[a][level].add(center)
                pointsA[a] = level
        for b in self.B:
            #if self.dist(center,b) < min(self.nearestA[b], self.nearestB[b]):
            if self.dist(center,b) < self.nearestA[b]:
                level = -inf
                if self.dist(center,b) > 0:
                    level = ceil(log(self.dist(center,b), self.base))
                    self.levels.add(level)
                bucketsB[level].add(b)
                self.lookup[b][level].add(center)
                pointsB[b] = level
        if len(bucketsA) > 0 and len(bucketsB) > 0:
            self.clusters[center] = Cluster(pointsA, pointsB,center, bucketsA, bucketsB, self.p, self.base, self.delta, self.diam, self)

    def compute_clusterB(self, center, i):
        #if VERBOSE: print('.', end='', flush=True)
        bucketsA, bucketsB = defaultdict(set), defaultdict(set)
        pointsA, pointsB = dict(), dict()
        for a in self.A:
            #if self.dist(a,center) < min(self.nearestB[a], self.nearestA[a]):
            if self.dist(a,center) < self.nearestB[a]:
                level = -inf
                if self.dist(a,center) > 0:
                    level = ceil(log(self.dist(a,center),self.base))
                    self.levels.add(level)
                bucketsA[level].add(a)
                self.lookup[a][level].add(center)
                pointsA[a] = level
        for b in self.B:
            #for b in self.layersB[i]:
            #if self.dist(center,b) < min(self.nearestB[b],self.nearestA[b]):
            #if (i == 1 and b in self.layersB[i]) or (i != 1 and self.dist(center,b) < self.nearestB[b]):
            if (i != 1 and self.dist(center,b) < self.nearestB[b]):
                level = -inf
                if self.dist(center,b) > 0:
                    level = ceil(log(self.dist(center,b), self.base))
                    self.levels.add(level)
                bucketsB[level].add(b)
                self.lookup[b][level].add(center)
                pointsB[b] = level
        if len(bucketsA) > 0 and len(bucketsB) > 0:
            self.clusters[center] = Cluster(pointsA, pointsB,center, bucketsA, bucketsB, self.p, self.base, self.delta, self.diam, self)

    def init_min_slack_heap(self):
        self.slack_heap = MinHeap()
        for center, cluster in self.clusters.items():
            cluster.init_max_heaps()
            self.slack_heap.insert(cluster, self._max_slack)
        for a in self.A: a.inserted = True
        for b in self.B: b.inserted = True

    def bcp(self):
        cluster = self.slack_heap.findmin()
        return cluster.min_slack_edge

    def activate_all_A(self, exclude=set()):
        for a in self.inactivated_A:
            if a not in exclude:
                for level in self.lookup[a]:
                    for center in self.lookup[a][level]:
                        if center not in self.clusters: continue 
                        cluster = self.clusters[center]
                        cluster.updateA(a)
        self.inactivated_A = set(exclude)

    def removeA(self, a, lazy=False):
        a.inserted = False
        self.inactivated_A.add(a)
        for level in self.lookup[a]:
            for center in self.lookup[a][level]:
                if center not in self.clusters: continue
                cluster = self.clusters[center]
                cluster.heapA.remove(a)
                if not lazy and cluster.min_slack_edge.a == a:
                    edge = cluster.find_BCP()
                    self.slack_heap.changepriority(cluster, edge.slack)

    def removeB(self, b, lazy=False):
        b.inserted = False
        if b in self.activated_B:
            self.activated_B.remove(b)
            for level in self.lookup[b]:
                for center in self.lookup[b][level]:
                    if center not in self.clusters: continue
                    cluster = self.clusters[center]
                    cluster.heapB.remove(b)
                    if not lazy and cluster.min_slack_edge.b == b:
                        edge = cluster.find_BCP()
                        self.slack_heap.changepriority(cluster, edge.slack)

    def activate_only_B(self, B):
        rmv_pts = self.activated_B - B
        for b in rmv_pts: self.removeB(b, lazy=True)
        self.activated_B = set(B)
        for b in B:
            b.inserted = True
            self.updateB(b, lazy=True)
        for center, cluster in self.clusters.items():
            edge = cluster.find_BCP()
            self.slack_heap.changepriority(cluster, edge.slack)

    def find_NN(self, b):
        best_edge = Edge()
        if b in self.clusters:
            cluster = self.clusters[b]
            best_edge = cluster.find_NN()
        else:
            for level in self.lookup[b]:
                for center in self.lookup[b][level]:
                    if center not in self.clusters: continue
                    cluster = self.clusters[center]
                    edge = cluster.find_NN_b(b)
                    if edge.slack <= best_edge.slack:
                        best_edge = edge
        return best_edge

    def updateB(self, b, lazy=False):
        self.activated_B.add(b)
        b.inserted = True
        for level in self.lookup[b]:
            for center in self.lookup[b][level]:
                if center not in self.clusters: continue 
                cluster = self.clusters[center]
                cluster.updateB(b)
                if not lazy:
                    edge = cluster.find_BCP()
                    self.slack_heap.changepriority(cluster, edge.slack)

    def compute_cluster_dist(self):
        self.distC = dict()
        self.proxyDistC = dict()
        for a in self.A:
            for b in self.B:
                self.distC[(a,b)] = self.diam
                self.proxyDistC[(a,b)] = ceil(self.diam/self.delta)
        for center, cluster in self.clusters.items():
            for a in cluster.A:
                for b in cluster.B:
                    self.distC[(a,b)] = min(cluster.distC(a,b), self.distC[(a,b)])
                    self.proxyDistC[(a,b)] = min(cluster.proxyDistC(a,b), self.proxyDistC[(a,b)])

    def print_cluster_heaps(self, visitedA, visitedB):
        min_slack = ceil(pow(self.diam, self.p)/self.delta)
        for center, cluster in self.clusters.items():
            cluster.print_heap_contents()
            for a in cluster.A:
                for b in cluster.B:
                    if a not in visitedA and b in visitedB:
                        slack = cluster.proxyDistC(a,b) - a.dual_weight - b.dual_weight + b.len_path
                        min_slack = min(slack, min_slack)


    def print_clusters(self):
        print("--------------------------------")
        for cluster in self.clusters:
            self.clusters[cluster].print_cluster()
        print("--------------------------------")

class Cluster:
    def __init__(self, A, B, center, bucketsA, bucketsB, p, base, delta, diam, parent):
        self.parent = parent
        self.A = A
        self.B = B
        self.center = center
        self.id = center.id
        self.bucketsA = bucketsA
        self.bucketsB = bucketsB
        self.p = p
        self.base = base
        self.delta = delta
        self.diam = diam
        self._max_slack = ceil(pow(self.diam, self.p)/self.delta)
        all_levels_A, all_levels_B = set(bucketsA.keys()), set(bucketsB.keys())
        self.all_levels = list(all_levels_A.union(all_levels_B))
        self.all_levels.sort()
        self.min_slack_edge = Edge(None, None, self._max_slack, self.center, None)


    def find_BCP(self):
        max_a, max_b = None, None
        min_slack = self._max_slack
        min_slack_edge = Edge(None, None, min_slack, self.center, None)
        for level in self.all_levels:
            level_a = self.heapA.max_item_cache.get(level)
            level_b = self.heapB.max_item_cache.get(level)
            if self.center.inserted:
                if self.center in self.A and level_b is not None:
                    slack = self.parent.cache_proxy_dist[level] - self.center.dual_weight - (level_b.dual_weight - level_b.len_path) + 1
                    if min_slack > slack:
                        min_slack = slack
                        min_slack_edge = Edge(self.center, level_b, slack, self.center, level)
                if self.center in self.B and level_a is not None:
                    slack = self.parent.cache_proxy_dist[level] - level_a.dual_weight - (self.center.dual_weight - self.center.len_path) + 1
                    if min_slack > slack:
                        min_slack = slack
                        min_slack_edge = Edge(level_a, self.center, slack, self.center, level)
            if max_a is None or level_a is not None and level_a.dual_weight > max_a.dual_weight:
                max_a = level_a
            if max_b is None or level_b is not None and level_b.dual_weight - level_b.len_path > max_b.dual_weight - max_b.len_path:
                max_b = level_b
            if max_a is not None and max_b is not None:
                slack = self.parent.cache_proxy_diam[level] - max_a.dual_weight - (max_b.dual_weight - max_b.len_path) + 1

                if min_slack > slack:
                    min_slack = slack
                    min_slack_edge = Edge(max_a, max_b, slack, self.center, level)
        self.min_slack_edge = min_slack_edge
        self.min_slack = min_slack
        return min_slack_edge

    def find_NN(self):
        edge = Edge(None, self.center, self._max_slack, self.center, None)
        for level in self.all_levels:
            #print("level", level)
            a = self.heapA.max_item_cache.get(level)
            if edge.a is None: #or (a is not None and a.dual_weight > edge.a.dual_weight):
                edge.a = a
            if a is not None:
                slack = self.parent.cache_proxy_dist[level] - a.dual_weight - self.center.dual_weight + 1
                if slack < edge.slack:
                    edge.slack = slack
                    edge.a = a
                    #print("in find_NN", edge)
                    #print("w_a", edge.a.dual_weight, "w_b", edge.b.dual_weight)
                    #print("level", level)
                    #print("proxy dist", self.parent.cache_proxy_dist[level])
        return edge

    def find_NN_b(self, b):
        edge = Edge(None, b, self._max_slack, self.center, None)
        b_level = self.B[b]
        for level in self.all_levels:
            #print("level", level)
            a = self.heapA.max_item_cache.get(level)
            if edge.a is None: #or (a is not None and a.dual_weight > edge.a.dual_weight):
                edge.a = a
            if a is not None:
                if a == self.center:
                    slack = self.parent.cache_proxy_dist[max(level,b_level)] - a.dual_weight - b.dual_weight + 1
                    #print("level", level)
                    #print("proxy dist", self.parent.cache_proxy_dist[level])
                else:
                    slack = self.parent.cache_proxy_diam[max(level,b_level)] - a.dual_weight - b.dual_weight + 1
                    #print("level", level)
                    #print("proxy diam", self.parent.cache_proxy_diam[level])

                if slack < edge.slack:
                    edge.slack = min(edge.slack, slack)
                    edge.a = a
                    #print("in find_NN_b", edge)
                #print("w_a", edge.a.dual_weight, "w_b", edge.b.dual_weight)
        return edge

    def insertB(self, b):
        level = self.B[b]
        self.heapB.insert(b, b.dual_weight - b.len_path, level)
                
    #def removeA(self, a):
        #self.heapA.remove(a)
        #if self.min_slack_edge.a == a:
            #self.find_BCP()

    #def removeB(self, b):
        #self.heapB.remove(b)
        #if self.min_slack_edge.b == b:
            #self.find_BCP()

    def updateB(self, b):
        b.inserted = True
        if b in self.heapB:
            self.heapB.changepriority(b, b.dual_weight - b.len_path)
        elif b in self.B:
            level = self.B[b]
            self.heapB.insert(b, b.dual_weight - b.len_path, level)
        # TODO: check if b is part of the new best edge and update if needed

    def updateA(self, a):
        a.inserted = True
        if a in self.heapA:
            self.heapA.changepriority(a, a.dual_weight)
        elif a in self.A:
            level = self.A[a]
            self.heapA.insert(a, a.dual_weight, level)
        # TODO: check if a is part of the new best edge and update if needed

    def init_max_heaps(self):
        self.heapA, self.heapB = MaxLevelHeap(), MaxLevelHeap()
        for level, bucket in self.bucketsA.items():
            for a in bucket:
                #self.heapA.insert(a, a.dual_weight, level)
                self.heapA.insert(a, 0, level)
        for level, bucket in self.bucketsB.items():
            for b in bucket:
                #self.heapB.insert(b, b.dual_weight - b.len_path,level)
                self.heapB.insert(b, 0, level)

    '''
    Computes a sorted list of the non-empty levels
    '''
    def compute_active_levels(self):
        self.active_levels_A = self.heapA.levels
        self.active_levels_B = self.heapB.levels
        active_levels = self.active_levels_A.union(self.active_levels_B)
        self.active_levels = list(levels)
        self.active_levels.sort()

    '''
    Cluster distance within this cluster
    '''
    def distC(self, ptA, ptB):
        level = max(self.A[ptA], self.B[ptB])
        if level is -inf: return 0
        #return 2*pow(self.base, level)
        if ptA.id == self.id or ptB.id == self.id:
            return self.parent.cache_level_dist[level]
        else:
            return 2*self.parent.cache_level_dist[level]

    def proxyDistC(self, ptA, ptB):
        level = max(self.A[ptA], self.B[ptB])
        if level is -inf: return 0
        if ptA.id == self.id or ptB.id == self.id:
            return self.parent.cache_proxy_dist[level]
        else:
            return self.parent.cache_proxy_diam[level]

    def print_cluster(self):
        print("Cluster", self.center.id, end=": ")
        print(" ".join([str(a.id) for a in self.A]), end=" + ")
        print(" ".join([str(b.id) for b in self.B]))

    def print_heap_contents(self):
        #print("in print_heap_contents() max val = ", self.heapB.max_val)
        if not self.heapA.isEmpty() and not self.heapB.isEmpty():
            print("Cluster", self.center.id, ":")
            print("Heap A:", self.heapA)
            print("Heap B:", self.heapB)
            #print("in print_heap_contents() max val = ", self.heapB.max_val)

