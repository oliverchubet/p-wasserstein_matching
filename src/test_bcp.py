import unittest
from decomp import Decomposition
from collections import defaultdict
from math import log, ceil
import ot
from constants import *
from distribution import generate_points
import math
import utils
import csv
from point import Point

distance_function = utils.dist

class TestClustering(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.p = 2
        cls.delta = 0.01
        cls.k = 2
        cls.n = 30
        cls.dim = 2
        cls.base = 1.01
        cls.distribution = "Normal"
        cls.A, cls.B, cls.masses_A, cls.masses_B = generate_points(cls.n,cls.dim,cls.distribution)
        cls.decomp = Decomposition(cls.A, cls.B, distance_function, cls.p, cls.delta, cls.k, cls.base, cls.dim)
        cls.clusters = cls.decomp.clusters
        cls.decomp.compute_cluster_dist()
        cls.distC = cls.decomp.distC
        cls.diam = cls.decomp.diam
        cls.proxyDistC = cls.decomp.proxyDistC

    def test_proxy_dist_in_cluster(self):
        for cluster in self.clusters.values():
            for a in cluster.A:
                for b in cluster.B:
                    pd = cluster.proxyDistC(a,b)
                    self.assertEqual(pd, ceil(pow(cluster.distC(a,b), self.p)/self.delta))



    def test_init_min_slack_heap(self):
        for b in self.B:
            b.len_path = 0
        self.decomp.init_min_slack_heap(self.A, self.B)
        edge = self.decomp.weighted_BCP()
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x])
        self.assertEqual(self.proxyDistC[correct_edge], self.proxyDistC[(edge.a, edge.b)])

    def test_weighted_BCP_after_updates(self):
        for b in self.B:
            b.len_path = 0
        self.decomp.init_min_slack_heap(self.A, self.B)
        edge = self.decomp.weighted_BCP()
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x])
        self.assertEqual(self.proxyDistC[correct_edge], self.proxyDistC[(edge.a, edge.b)])
        b = self.B[0]
        b.len_path += 10
        self.decomp.updateB(b)
        edge = self.decomp.weighted_BCP()
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x] + x[1].len_path)
        if self.proxyDistC[correct_edge]+correct_edge[1].len_path != self.proxyDistC[(edge.a, edge.b)]+edge.b.len_path:
            print("correct:", "".join(["(",str(correct_edge[0].id),",", str(correct_edge[1].id),");"]),"slack =", self.proxyDistC[correct_edge]+correct_edge[1].len_path)
            print("returned:", edge)
            print(self.decomp.slack_heap)
        self.assertEqual(self.proxyDistC[correct_edge]+correct_edge[1].len_path, self.proxyDistC[(edge.a, edge.b)]+edge.b.len_path)

    def test_weighted_BCP_after_updates_2(self):
        for b in self.B:
            b.len_path = 0
        self.decomp.init_min_slack_heap(self.A, self.B)
        edge = self.decomp.weighted_BCP()
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x])
        self.assertEqual(self.proxyDistC[correct_edge], self.proxyDistC[(edge.a, edge.b)])
        b = self.B[0]
        b.len_path += 10
        self.decomp.updateB(b)
        edge = self.decomp.weighted_BCP()
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x] + x[1].len_path)
        self.assertEqual(self.proxyDistC[correct_edge]+correct_edge[1].len_path, self.proxyDistC[(edge.a, edge.b)]+edge.b.len_path)
        edge.b.len_path += 10
        self.decomp.updateB(edge.b)
        edge = self.decomp.weighted_BCP()
        #for center, cluster in self.clusters.items():
        #    print("id:", center.id,"min slack:", cluster.min_slack)
        #print(self.decomp.slack_heap)
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x] + x[1].len_path)
        if self.proxyDistC[correct_edge]+correct_edge[1].len_path != self.proxyDistC[(edge.a, edge.b)]+edge.b.len_path:
            print("correct:", "".join(["(",str(correct_edge[0].id),",", str(correct_edge[1].id),");"]),"slack =", self.proxyDistC[correct_edge]+correct_edge[1].len_path)
            print("returned:", edge)
            print(self.decomp.slack_heap)
        self.assertEqual(self.proxyDistC[correct_edge]+correct_edge[1].len_path, self.proxyDistC[(edge.a, edge.b)]+edge.b.len_path)

    def test_weighted_BCP_per_cluster(self):
        for b in self.B:
            b.len_path = 0
        self.decomp.init_min_slack_heap(self.A, self.B)
        for cluster in self.clusters.values():
            #for a in cluster.A:
            #    for b in cluster.B:
            #        print(a.id, b.id, ":", cluster.proxyDistC(a,b) + b.len_path)
            edge = cluster.weighted_BCP()
            #print("returned", edge.a.id, edge.b.id)
            correct_edge = min([(a,b) for a in cluster.A for b in cluster.B ], key=lambda x: cluster.proxyDistC(*x) + x[1].len_path)
            #print("correct", correct_edge[0].id, correct_edge[1].id)
            self.assertEqual(cluster.proxyDistC(*correct_edge)+correct_edge[1].len_path, cluster.proxyDistC(edge.a, edge.b)+edge.b.len_path)

    def test_weighted_BCP_per_cluster_after_updates(self):
        for b in self.B:
            b.len_path = 0
        self.decomp.init_min_slack_heap(self.A, self.B)
        edge = self.decomp.weighted_BCP()
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x])
        self.assertEqual(self.proxyDistC[correct_edge], self.proxyDistC[(edge.a, edge.b)])
        b = self.B[0]
        b.len_path += 10
        self.decomp.updateB(b)
        #print("First b should have len_path 10")
        #for cluster in self.clusters.values():
        #    cluster.print_heap_contents()
        edge = self.decomp.weighted_BCP()
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x] + x[1].len_path)
        self.assertEqual(self.proxyDistC[correct_edge]+correct_edge[1].len_path, self.proxyDistC[(edge.a, edge.b)]+edge.b.len_path)

    def test_weighted_BCP_per_cluster_after_updates_2(self):
        for b in self.B:
            b.len_path = 0
        self.decomp.init_min_slack_heap(self.A, self.B)
        edge = self.decomp.weighted_BCP()
        correct_edge = min([(a,b) for a in self.A for b in self.B ], key=lambda x: self.proxyDistC[x] + x[1].len_path)
        self.assertEqual(self.proxyDistC[correct_edge]+correct_edge[1].len_path, self.proxyDistC[(edge.a, edge.b)]+edge.b.len_path)
        edge.b.len_path += 10
        self.decomp.updateB(edge.b)
        #print("Two b's should have len_path 10")
        for cluster in self.clusters.values():
            #cluster.print_heap_contents()
        #    for a in cluster.A:
        #        for b in cluster.B:
        #            print(a.id, b.id, ":", cluster.proxyDistC(a,b) + b.len_path)
            #print("returned", edge.a.id, edge.b.id)
            correct_edge = min([(a,b) for a in cluster.A for b in cluster.B ], key=lambda x: cluster.proxyDistC(*x) + x[1].len_path)
            edge = cluster.weighted_BCP()
            #print("returned", edge.a.id, edge.b.id, cluster.proxyDistC(a,b)+edge.b.len_path)
            if cluster.proxyDistC(*correct_edge)+correct_edge[1].len_path != cluster.proxyDistC(edge.a, edge.b)+edge.b.len_path:
                print("correct:", "".join(["(",str(correct_edge[0].id),",", str(correct_edge[1].id),");"]),"slack =", cluster.proxyDistC(*correct_edge)+correct_edge[1].len_path)
                print("returned:", edge)
                cluster.print_heap_contents()
            self.assertEqual(cluster.proxyDistC(*correct_edge)+correct_edge[1].len_path, cluster.proxyDistC(edge.a, edge.b)+edge.b.len_path)

    def test_remove(self):
        for b in self.B:
            b.len_path = 0
        self.decomp.init_min_slack_heap(self.A, self.B)
        edge = self.decomp.weighted_BCP()
        self.decomp.removeA(edge.a)
        edge2 = self.decomp.weighted_BCP()
        self.assertNotEqual(edge.a, edge2.a)
