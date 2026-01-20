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
        #if VERBOSE: print("Compute decomposition")
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

        self.compute_samples()
        for i in range(self.k-1, 0, -1):
            ##if VERBOSEPLUS: print("Computing layer", i, " clusters", end=" ")
            self.compute_clusters(i)
            self.compute_nearest(i)
        ##if VERBOSEPLUS: print("\nComputing layer 0 clusters", end=" ")
        self.compute_clusters(0)
        ##if VERBOSEPLUS: print()
        self.cache_level_dist = {level: pow(base, level) for level in self.levels}
        self.cache_proxy_dist = {level: ceil(pow(self.cache_level_dist[level],p)/self.delta) for level in self.levels}
        self.cache_proxy_diam = {level: ceil(pow(2*self.cache_level_dist[level],p)/self.delta) for level in self.levels}

    '''
    Each point in A and B is assigned to a layer
    '''
    def compute_samples(self):
        prob = pow(self.n, -1.0/self.k)
        #for pt in chain(self.A, self.B):
        for pt in chain(self.A):
            for i in range(self.k):
                if i == self.k-1 or random.random() > prob:
                    self.layersA[i].add(pt)
                    break
        for pt in chain(self.B):
            for i in range(self.k):
                if i == self.k-1 or random.random() > prob:
                    self.layersB[i].add(pt)
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
            self.compute_clusterA(c)
        for c in self.layersB[i]:
            self.compute_clusterB(c)

    '''
    Adds points to the cluster if the center is closer than the nearest neighbor distance to the sample.
    Points are bucketed by levels determined by distance to the center.
    '''
    def compute_clusterA(self, center):
        #if VERBOSE: print('.', end='', flush=True)
        bucketsA, bucketsB = defaultdict(set), defaultdict(set)
        pointsA, pointsB = dict(), dict()
        for a in self.A:
            if self.dist(a,center) < self.nearestA[a]:
                level = -inf
                if self.dist(a,center) > 0:
                    level = ceil(log(self.dist(a,center),self.base))
                    self.levels.add(level)
                bucketsA[level].add(a)
                self.lookup[a][level].add(center)
                pointsA[a] = level
        for b in self.B:
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

    def compute_clusterB(self, center):
        #if VERBOSE: print('.', end='', flush=True)
        bucketsA, bucketsB = defaultdict(set), defaultdict(set)
        pointsA, pointsB = dict(), dict()
        for a in self.A:
            if self.dist(a,center) < self.nearestB[a]:
                level = -inf
                if self.dist(a,center) > 0:
                    level = ceil(log(self.dist(a,center),self.base))
                    self.levels.add(level)
                bucketsA[level].add(a)
                self.lookup[a][level].add(center)
                pointsA[a] = level
        for b in self.B:
            if self.dist(center,b) < self.nearestB[b]:
                level = -inf
                if self.dist(center,b) > 0:
                    level = ceil(log(self.dist(center,b), self.base))
                    self.levels.add(level)
                bucketsB[level].add(b)
                self.lookup[b][level].add(center)
                pointsB[b] = level
        if len(bucketsA) > 0 and len(bucketsB) > 0:
            self.clusters[center] = Cluster(pointsA, pointsB,center, bucketsA, bucketsB, self.p, self.base, self.delta, self.diam, self)

    '''
    Take the min slack edge from each cluster and put it into a min heap
    '''
    def init_min_slack_heap(self, A, B): #, tlevel): #threshhold=None):
        self.slack_heap = MinHeap()
        for center, cluster in self.clusters.items():
            cluster.init_max_heaps(A, B) #, tlevel)
            #cluster.print_heap_contents()
            edge = cluster.weighted_BCP()
            self.slack_heap.insert(cluster, edge.slack)
        for a in self.A:
            a.inserted = a in A
        for b in self.B:
            b.inserted = b in B

    def update_slack_heap(self):
        ##if VERBOSEPLUS: print("Updating min slack heap")
        for center, cluster in self.clusters.items():
            edge = cluster.weighted_BCP()
            #if cluster.min_slack > edge.slack:
            self.slack_heap.changepriority(cluster, edge.slack)


    '''
    Returns the min slack edge but does not remove it from the heap
    '''
    def weighted_BCP(self):
        if TRACE: print("decomp.weighted_BCP")
        self.update_slack_heap()
        #print(self.slack_heap)
        if len(self.slack_heap) > 0:
            cluster = self.slack_heap.findmin()
            return cluster.min_slack_edge

    #def weighted_NN(self, x):
    #    min_slack_edge = None
    #    min_slack = self.diam
    #    for level in self.lookup[x]:
    #        for center in self.lookup[x][level]:
    #            if center not in self.clusters:
    #                continue
    #            edge = self.clusters[center].weighted_NN(x)
    #            if edge is None:
    #                continue
    #            if min_slack_edge is None or edge.slack < min_slack:
    #                min_slack_edge = edge
    #                min_slack = edge.slack
    #    return min_slack_edge

    '''
    Each loop of Djikstra's, some points will have their path length updated.
    Tell the clusters affected to update their heaps and min slack edge
    Updating the min slack heap
    '''
    def removeA(self, a): #, tlevel): #threshhold=None):
        #self.clusters[a.id].inserted = False
        a.inserted = False
        for level in self.lookup[a]:
            #if level <= tlevel:
            for center in self.lookup[a][level]:
                if center not in self.clusters:
                    continue
                cluster = self.clusters[center]
                cluster.removeA(a)
                #if cluster not in self.slack_heap:
                #    continue
                edge = cluster.weighted_BCP()
                #if edge is None:
                #    if cluster in self.slack_heap:
                #        self.slack_heap.remove(cluster)
                #else:
                self.slack_heap.changepriority(cluster, edge.slack)

    def removeB(self, b): #, tlevel): #threshhold=None):
        #self.clusters[b.id].inserted = False
        b.inserted = False
        for level in self.lookup[b]:
            #if level <= tlevel:
            for center in self.lookup[b][level]:
                if center not in self.clusters:
                    continue
                cluster = self.clusters[center]
                #if cluster not in self.slack_heap:
                #    continue
                cluster.removeB(b)
                edge = cluster.weighted_BCP() #threshhold=threshhold)
                #if edge is None:
                #    self.slack_heap.remove(cluster)
                #else:
                self.slack_heap.changepriority(cluster, edge.slack)

    def updateB(self, b): #, tlevel): #threshhold=None):
        for level in self.lookup[b]:
            #if level <= tlevel:
            for center in self.lookup[b][level]:
                if center not in self.clusters:
                    continue 
                cluster = self.clusters[center]
                cluster.updateB(b)
        self.update_slack_heap()

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
            #if DEBUG: print("Cluster", center.id)
            cluster.print_heap_contents()
            for a in cluster.A:
                for b in cluster.B:
                    if a not in visitedA and b in visitedB:
                        slack = cluster.proxyDistC(a,b) - a.dual_weight - b.dual_weight + b.len_path
                        #if DEBUG: print("a:", a.id,"b:", b.id, "proxyDist:", cluster.proxyDistC(a,b), "y(a):", a.dual_weight, "y(b):", b.dual_weight, "l(b):", b.len_path, "slack:", slack)
                        min_slack = min(slack, min_slack)
        #if DEBUG: print("min slack:", min_slack)


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
        self._proxy_diam = ceil(pow(self.diam, self.p)/self.delta)
        # caching the value 2*pow(base, level) to avoid function call overhead
        levelsA, levelsB = set(bucketsA.keys()), set(bucketsB.keys())
        self.levels = list(levelsA.union(levelsB))
        self.levels.sort()
        #self.threshold = self.levels[0]
        #self._cache_level_dist = {level: pow(base, level) for level in levels}
        #self._cache_proxy_dist = {level: ceil(pow(self._cache_level_dist[level],p)/self.delta) for level in levels}
        #self._cache_proxy_diam = {level: ceil(pow(2*self._cache_level_dist[level],p)/self.delta) for level in levels}
        #self._cache_proxy_dist = {level: pow(self._cache_level_dist[level], p) for level in levels}
        #self._cache_proxy_diam = {level: pow(2*self._cache_level_dist[level], p) for level in levels}
        self.inserted = False


    '''
    Once the max heaps and max weights are initialized, the min slack edge is computed
    '''

    def weighted_BCP(self):
        #if TRACE: print("cluster.weighted_BCP")
        ##if VERBOSEPLUS: print("\tLooking for BCP of cluster", self.center.id, "...")
        max_a, max_b = None, None
        min_slack = self._proxy_diam
        min_slack_edge = Edge(None, None, min_slack, self.center, None)
        #print("Cluster", self.center.id, end=" ")
        for level in self.levels:
            #if level < threshold:
            if DEBUGPLUS: print("\t\t\tlevel", level, end=" ")
            #level_a = self.heapA.findmax(level)
            #level_b = self.heapB.findmax(level)
            level_a = self.heapA.max_item_cache.get(level)
            level_b = self.heapB.max_item_cache.get(level)
            if DEBUGPLUS and max_a: print("max_a", max_a.id, end=" ")
            if DEBUGPLUS and max_b: print("max_b", max_b.id, end="") 
            if self.center.inserted:
                if self.center in self.A and level_b is not None:
                    slack = self.parent.cache_proxy_dist[level] - self.center.dual_weight - (level_b.dual_weight - level_b.len_path) + 1
                    if min_slack > slack:
                        min_slack = slack
                        min_slack_edge = Edge(self.center, level_b, slack, self.center, level)
                        #print(min_slack_edge)
                if self.center in self.B and level_a is not None:
                    slack = self.parent.cache_proxy_dist[level] - level_a.dual_weight - (self.center.dual_weight - self.center.len_path) + 1
                    if min_slack > slack:
                        min_slack = slack
                        min_slack_edge = Edge(level_a, self.center, slack, self.center, level)
            if max_a is None or level_a is not None and level_a.dual_weight > max_a.dual_weight:
                max_a = level_a
            if max_b is None or level_b is not None and level_b.dual_weight - level_b.len_path > max_b.dual_weight - max_b.len_path:
                max_b = level_b
            if DEBUGPLUS: print()
            if max_a is not None and max_b is not None:
                #slack = self.proxyDistC(max_a, max_b) - max_a.dual_weight - (max_b.dual_weight - max_b.len_path) + 1
                #if max_a.id == self.id or max_b.id == self.id:
                #    slack = self._cache_proxy_level_dist[level] - max_a.dual_weight - (max_b.dual_weight - max_b.len_path) + 1
                #else:
                slack = self.parent.cache_proxy_diam[level] - max_a.dual_weight - (max_b.dual_weight - max_b.len_path) + 1

                if DEBUGPLUS: print("\t\tFound", max_a.id, max_b.id, "slack:", slack)
                #print("slack =", slack, end="")
                if min_slack > slack:
                    min_slack = slack
                    min_slack_edge = Edge(max_a, max_b, slack, self.center, level)
        self.min_slack_edge = min_slack_edge
        self.min_slack = min_slack
        #if DEBUGPLUS: print("min slack edge:", self.min_slack_edge)
        #print("center", self.id, "min slack edge:", self.min_slack_edge)
        if DEBUGPLUS: print()
        #print("Cluster", self.id, self.min_slack_edge)
        return min_slack_edge

#    def weighted_NN(self,b):
#        max_a = None
#        min_slack = self.diam
#        min_slack_edge = None
#        for level in self.levels:
#            #max_a = self.heapA.findmax(level)
#            max_a = self.heapA.max_item_cache.get(level)
#            if max_a is not None:
#                #slack = self.proxyDistC(max_a, b) - max_a.dual_weight - (b.dual_weight - b.len_path)
#                if b.id == self.id or max_a.id == self.id:
#                    slack = self._cache_proxy_level_dist[level] - max_a.dual_weight - (b.dual_weight - b.len_path)
#                else:
#                    slack = 2*self._cache_proxy_level_dist[level] - max_a.dual_weight - (b.dual_weight - b.len_path)
#                if min_slack > slack:
#                    min_slack = slack
#                    min_slack_edge = Edge(max_a, b, slack, self.center, level)
#        return min_slack_edge

    def removeA(self, a):
        #level = self.A[a]
        #if level < threshold:
        self.heapA.remove(a)

    def removeB(self, b):
        #level = self.B[b]
        #if level < threshold:
        self.heapB.remove(b)

    def insertB(self, b):
        level = self.B[b]
        #if level < threshold:
        self.heapB.insert(b, b.dual_weight - b.len_path, level)
                
    def updateB(self, b):
        ##if VERBOSEPLUS: print("Cluster", self.center.id, "updating", b.id)
        b.inserted = True
        level = self.B[b]
        self.levelsB.add(level)
        self.compute_levels()
        #if level < threshold:
        if b in self.heapB:
            self.heapB.changepriority(b, b.dual_weight - b.len_path)
        else:
            self.heapB.insert(b, b.dual_weight - b.len_path, level)
        #self.weighted_BCP()
        #self.print_heap_contents()

    '''
    Insert points into heaps by level keyed by their weight
    Initializes self.heapsA and self.heapsB (each level maps to a max heap keyed by weights)
    '''
    def init_max_heaps(self, A, B): #, tlevel):
        self.heapA, self.heapB = MaxLevelHeap(), MaxLevelHeap()
        for level, bucket in self.bucketsA.items():
            #if level <= tlevel:
            for a in bucket:
                if a in A:
                    self.heapA.insert(a, a.dual_weight, level)
        for level, bucket in self.bucketsB.items():
            #if level < tlevel:
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
        #if None in levels: levels.remove(None)
        self.levels = list(levels)
        self.levels.sort()
        #self.levels = [None] + levels

    #def compute_levels(self):
    #    levelsA

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
