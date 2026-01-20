from constants import *
from distribution import generate_points
from utils import dist
from decomp import Decomposition
from itertools import chain
from tree import Tree
import copy
from math import ceil, log, inf
import ot

#THRESHOLD = True

class Wasserstein:
    def __init__(self, A, B, dist, p=1, delta=0.01, base=1.01, dim=2, k=2):
        self.A = A
        self.B = B
        self.n = len(A)
        self.dist = dist
        self.p = p
        self.diam = 2*pow(2, dim)
        self.min_cost = 0
        self.cost_using_clustering = 0
        self.decomp = Decomposition(A, B, dist, p, delta, k=k, base=base, dim=dim)
        self.freeB = set([b for b in self.B])
        self.visitedA = set()
        self.visitedB = set([b for b in self.freeB])
        self.unvisitedA = set([a for a in self.A])
        self.matched = 0
        self.phases = 0
        self.edges_seen = 0
        self.delta = delta
        self.base = base
        #self.tlevel = ceil(log(delta, self.base))
        #if THRESHOLD is False:
            #self.tlevel = ceil(log(self.diam, self.base))

    '''
    Computes min cost matching using p-th power distances
    '''
    def compute_pWasserstein(self):
        #if VERBOSE: print("\nComputing pWasserstein distance...")
        count, len_path = 0, 0
        while self.matched < self.n and count < self.n: 
            count += 1
            #if VERBOSE: print("###################################")
            print("Iteration", count,":")
            self.run_phase()
            if DEBUG: self.print_slack_matrix()
            if DEBUG: self.check_all_feasible()
        assert(self.is_matching())
        self.check_all_feasible()
        self.compute_cost()

    def run_phase(self):
        len_path = self.djikstra()
        self.update_dual_weights(len_path)
        self.partial_dfs()

    def compute_cost(self):
        cluster_cost = sum([pow(self.distC[(a,a.match)],self.p) for a in self.A])
        proxy_cost = sum([self.proxyDistC[(a,a.match)] for a in self.A])
        self.cost_using_clustering = pow(cluster_cost, 1.0/self.p)
        self.cost_using_proxy = pow(proxy_cost, 1.0/self.p)
        cost = sum([pow(dist(a,a.match),self.p) for a in self.A])
        self.min_cost = pow(cost, 1.0/self.p)
        num_preserved = 0
        for a in self.A:
            #print(a.id, "--", a.match.id, ":", end=" ")
            #print("dist =", dist(a,a.match),";", end=" ")
            #print("distC =", self.distC[(a,a.match)], end="; ")
            ratio = self.distC[(a,a.match)]/(dist(a,a.match))
            #print("ratio =", ratio)
            if ratio <= self.base: num_preserved += 1
        print("num preserved:", num_preserved)

    '''
    Uses the shortest path length returned by Djikstra's to update the dual weights
    '''
    def update_dual_weights(self, shortest_path):
        if VERBOSE: print("\tUpdating dual weights...")
        #threshold = pow(self.base, self.tlevel)
        if shortest_path is not None:
            for b in self.visitedB:
                b.dual_weight += shortest_path - b.len_path
                #if self.base*self.base*b.dual_weight > threshold:
                    #self.tlevel = ceil(log(b.dual_weight,self.base)) + 2
            for a in self.visitedA:
                a.dual_weight += a.len_path - shortest_path
        ###if VERBOSEPLUS: self.print_dual_weights()

    def check_all_feasible(self):
        if VERBOSE: print("\tChecking feasibility...")
        self.compute_cluster_distance()
        feasible = True
        for a in self.A:
            for b in self.B:
                feasible = feasible and self.is_feasible(a,b)
        assert(feasible)

    def is_feasible(self, a, b):
        proxy_dist = self.proxyDistC[(a,b)]
        if a.dual_weight + b.dual_weight > proxy_dist + 1:
            print("Not feasible...")
            print("proxy dist =", proxy_dist)
            print("".join(["y(",str(a.id),") = ", str(a.dual_weight)]))
            print("".join(["y(",str(b.id),") = ", str(b.dual_weight)]))
            return False

        if a.match is b:
            if a.dual_weight + b.dual_weight != proxy_dist:
                print("Not admissible...")
                print("proxy dist = ", proxy_dist)
                print("".join(["y(",str(a.id),") = ", str(a.dual_weight)]))
                print("".join(["y(",str(b.id),") = ", str(b.dual_weight)]))
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
    def djikstra(self): #threshhold=None):
        if VERBOSE: print("\tRunning Djikstra's...")
        self.visitedA = set()
        self.visitedB = set([b for b in self.freeB])
        self.unvisistedA = set([a for a in self.A])
        for b in self.B:
            b.len_path = None
            b.inserted = False
        for b in self.freeB:
            b.len_path = 0
            b.inserted = True
        for a in self.A:
            a.len_path = None
            a.inserted = True
        self.init_shortest_path_tree()
        self.decomp.init_min_slack_heap(self.A, self.freeB) #, self.tlevel)
        shortest_path = None
        count = 0
        while len(self.decomp.slack_heap) > 0 and count < self.n * self.n:
            count += 1
            if DEBUGPLUS: print(self.decomp.slack_heap)
            if DEBUGPLUS: self.decomp.print_cluster_heaps(self.visitedA, self.visitedB)
            edge = self.decomp.weighted_BCP()
            a, b = edge.a, edge.b
            if DEBUGPLUS and a is not None: print("l(", a.id, ") =", edge.slack)
            a.len_path = edge.slack #+ b.len_path
            if VERBOSEPLUS: print("\t\tConnect",b.id,"-",a.id, "weight:", edge.slack)
            self.shortest_path_tree.connect(a, b)
            self.visitedA.add(a)
            self.visitedB.add(b)
            self.decomp.removeA(a) #, self.tlevel) #threshhold=threshhold)
            if a.match is None:
                aug_path = self.shortest_path_tree.path(a)
                shortest_path = a.len_path
                self.augment(aug_path)
                return shortest_path
            else:
                b = a.match
                b.len_path = a.len_path
                #if DEBUG: print("l(", b.id, ") =", a.len_path)
                if VERBOSEPLUS: print("\t\tConnect",a.id,"-",b.id)
                self.shortest_path_tree.connect(b, a)
                self.visitedB.add(b)
                #if VERBOSEPLUS: print("\t\tRemove", a.id)
                self.decomp.removeA(a) #, self.tlevel) #threshhold=threshhold)
                #if VERBOSEPLUS: print("\t\tUpdate", b.id)
                if VERBOSEPLUS: self.print_slack_matrix(partial=True)
                self.decomp.updateB(b) #, self.tlevel) #threshhold=threshhold)
        #assert(count >= 1)
        #print("Shortest path len", shortest_path)
        return shortest_path

    def partial_dfs(self):
        #print("partial_dfs")
        #if VERBOSE: print("\tRunning partial DFS...")
        self.decomp.init_min_slack_heap(self.unvisitedA, []) #, self.tlevel) #threshhold=threshhold)
        pts = [b for b in self.freeB]
        for b in pts:
            #print("\tb is", b.id)
            aug_path = ["root", b]
            cur_b = b
            while True:
                self.decomp.updateB(cur_b) #, self.tlevel)
                edge = self.decomp.weighted_BCP()
                if edge is not None and edge.slack < 1:
                    #print("\t\t", edge)
                    #a,b = edge.a, edge.b
                    aug_path.append(edge.a)
                    #aug_path.append(b)
                    self.decomp.removeA(edge.a) #, self.tlevel)
                    self.decomp.removeB(cur_b) #, self.tlevel)
                    #self.decomp.removeB(b)
                    if edge.a.match is None:
                        #print("aug_path", aug_path[0], [x.id for x in aug_path[1:]])
                        self.augment(reversed(aug_path))
                    else:
                        cur_b = a.match
                        aug_path.append(cur_b)
                else:
                    aug_path.pop()
                    aug_path.pop()
                    if len(aug_path) <= 1: break
                    cur_b = aug_path[-1]

        #while len(self.decomp.slack_heap) > 0:
        #    edge = self.decomp.weighted_BCP()
        #    if edge.slack > 0: return
        #    self.shortest_path_tree.connect(edge.a, edge.b)
        #    self.decomp.removeA(edge.a)
        #    if edge.a.match is None:
        #        aug_path = self.shortest_path_tree.path(edge.a)
        #        self.decomp.removeA(edge.a)
        #        self.augment(aug_path)
        #    else:
        #        b = edge.a.match
        #        self.shortest_path_tree.connect(b, edge.a)
        #        self.visitedB.add(b)
        #        self.decomp.updateB(b)

    '''
    Updates the matching given the path returned by Djikstra's
    '''
    def augment(self, path):
        print("x", end="", flush=True)
        #if VERBOSE: print("\tAugmenting path", "-".join([str(pt.id) for pt in path[:-1]]), "...")
        it = iter(path)
        a = next(it)
        b = next(it)
        while True: #a != "root":
            self.phases += 1
            #b = next(it)
            if a in self.unvisitedA: self.unvisitedA.remove(a)
            b.dual_weight -= 1
            self.decomp.removeB(b) #, self.tlevel)
            if a == "root": break
            a.match = b
            a = next(it)
            if a == "root": break
            b = next(it)
            #self.shortest_path_tree.remove(a)
            #self.shortest_path_tree.remove(b)
            #print("b", b.id)
        self.freeB.remove(b)
        self.matched += 1
        #if DEBUG: self.print_dual_weights()

    def compute_cluster_distance(self):
        #if VERBOSE: print("Computing cluster distances...")
        self.distC = dict()
        self.proxyDistC = dict()
        for a in self.A:
            for b in self.B:
                self.distC[(a,b)] = self.diam 
                self.proxyDistC[(a,b)] = ceil(self.diam/self.delta)
        for center, cluster in self.decomp.clusters.items():
            for a in cluster.A:
                for b in cluster.B:
                    self.distC[(a,b)] = min(cluster.distC(a,b), self.distC[(a,b)])
                    self.proxyDistC[(a,b)] = min(cluster.proxyDistC(a,b), self.proxyDistC[(a,b)])
        #for a in self.A:
        #    for b in self.B:
        #        self.slackMatrix[(a,b)] = self.proxyDistC[(a,b)] - a.dual_weight - b.dual_weight

    def print_dual_weights(self):
        print("Dual weights:")
        print("+ A points:")
        for a in self.A:
            print("".join(["\ty(",str(a.id),") = ", str(a.dual_weight), "; l(",str(a.id),") = ",str(a.len_path)]))
        print("+ B points:")
        for b in self.B:
            print("".join(["\ty(",str(b.id),") = ", str(b.dual_weight), "; l(",str(b.id),") = ",str(b.len_path)]))

    def print_matching(self):
        print("Matched edges:")
        for a in self.A:
            if a.match is not None:
                print("".join(["\t(",str(a.id),", ",str(a.match.id),")"]))

    def print_slack_heap_edges(self):
        if self.decomp.slack_heap is not None and len(self.decomp.slack_heap) > 0:
            print("Cluster min slack edges:")
            for cluster in self.decomp.slack_heap._itemmap.keys():
                print("\tid:", cluster.center.id, "edge:", cluster.min_slack_edge.a.id, cluster.min_slack_edge.b.id, "slack:", cluster.min_slack_edge.slack)
        else:
            print("\tMin slack heap is empty...")

    def print_proxy_dist_matrix(self):
        self.compute_cluster_distance()
        for b in self.B:
            print("\t", b.id, end="")
        for a in self.A:
            print("\n", a.id, end="")
            for b in self.B:
                print("\t",self.proxyDistC[(a,b)],end="")
        print()

    def print_slack_matrix(self, partial=False):
        self.compute_cluster_distance()
        print("Slack Matrix")
        print("\t\t",end="")
        for b in self.B:
            print("\t", b.id, end="")
        print("\n\t\t y", end="")
        for b in self.B:
            print("\t",b.dual_weight, end="")
        print("\n\t y\t l", end="")
        for b in self.B:
            print("\t",b.len_path, end="")
        for a in self.A:
            print("\n", a.id, end="")
            print("\t",a.dual_weight, end="")
            print("\t",a.len_path, end="")
            for b in self.B:
                if not partial or a.len_path is None and b.len_path is not None:
                    print("\t",self.proxyDistC[(a,b)] - a.dual_weight - b.dual_weight + b.len_path,end="")
                else:
                    print("\t --",end="")
        print()
        print()

    def is_matching(self):
        for a in self.A:
            if a.match is None:
                print("Not matching?")
                self.print_matching()
                return False
        return True

    def print_matching(self):
        for a in self.A:
            if a.match:
                print(a.id, a.match.id)
            else:
                print(a.id, None)

if __name__ == "__main__":
    n = 100
    base = 1.1
    delta = 0.02
    p = 2
    A, B, masses_A, masses_B = generate_points(n,p,"Uniform")
    distance_function = dist
    wasserstein = Wasserstein(A, B, distance_function, p, delta=delta, base=base)
    wasserstein.compute_pWasserstein()
    dist_matrix = [ [ pow(dist(a,b),2) for a in A] for b in B ] 
    cluster_dist_matrix = [[ pow(wasserstein.distC[(a,b)],2) for a in A] for b in B]
    
    real_cost = pow(ot.emd2(masses_A, masses_B, dist_matrix), 1.0/2)
    cost = wasserstein.min_cost
    cluster_cost = wasserstein.cost_using_clustering
    cluster_emd_cost = pow(ot.emd2(masses_A, masses_B, cluster_dist_matrix), 1.0/2)
    cluster_ratio = cluster_cost / real_cost
    cluster_emd_ratio = cluster_emd_cost / real_cost
    ratio = cost / real_cost
    print()

    print("our matching cost: ", cost)
    print("\"real\" matching cost: ", real_cost)
    print("ratio: ", ratio)

    print()
    print("our cluster matching cost: ", cluster_cost)
    print("\"real\" matching cost: ", cluster_emd_cost)

    print()
    print("our cluster ratio: ", cluster_ratio)
    print("\"real\" cluster ratio: ", cluster_emd_ratio)


