import unittest
from skiplist import SkipList
from math import inf
from random import randint

SHAPES = True
DEBUG = False

class TestSkipList(unittest.TestCase):

    def test_skiplist_init(self):
        skeap = SkipList()
        node = skeap.find(-inf)
        self.assertEqual(node.key, -inf)
        node = skeap.find(0)
        self.assertEqual(node.key, -inf)
        node = skeap.find(1000)
        self.assertEqual(node.key, -inf)

    def test_skiplist_increase_height(self):
        skeap = SkipList()
        self.assertEqual(skeap.height, 2)
        skeap._increase_height()
        self.assertEqual(skeap.height, 3)
        skeap._increase_height()
        skeap._increase_height()
        skeap._increase_height()
        skeap._increase_height()
        skeap._increase_height()
        skeap._increase_height()
        skeap._increase_height()
        skeap._increase_height()
        self.assertEqual(skeap.height, 11)
        if SHAPES and DEBUG: skeap.print_shape()

    def test_skiplist_insert(self):
        skeap = SkipList()
        n = 50
        if SHAPES and DEBUG: 
            skeap.print_shape()
            n = 25
        for i in range(n):
            key, weight = randint(10,99), randint(10,99)
            if DEBUG: print(i, ";", key, weight, end=" ")
            skeap.insert(key, weight)
            if DEBUG:
                start = skeap.head
                cur = start
                for j in range(skeap.height):
                    while cur.next is not None:
                        print(cur.key, end=" ")
                        cur = cur.next
                    print()
                    start = start.down
                    cur = start
            if SHAPES and DEBUG: skeap.print_shape(weightues=True)
        if SHAPES and not DEBUG: skeap.print_shape()


    def test_skiplist_incremental_insert(self):
        skeap = SkipList()
        n = 50
        if SHAPES and DEBUG:
            skeap.print_shape()
            n = 25
        for i in range(10,10+n):
            key, weight = i, randint(10,99)
            skeap.insert(key, weight)
            if SHAPES and DEBUG: skeap.print_shape()
        if SHAPES and not DEBUG: skeap.print_shape()

    def test_skiplist_keys_in_order(self):
        skeap = SkipList()
        n = 500
        if SHAPES and DEBUG:
            skeap.print_shape()
            n = 25
        for i in range(10,10+n):
            key, weight = i, randint(10,99)
            skeap.insert(key, weight)
        start = skeap.head
        while start is not None:
            left, right = start, start.next
            while right is not None:
                self.assertTrue(left.key < right.key)
                left, right = right, right.next
            start = start.down

    def test_skiplist_remove(self):
        skiplist = SkipList()
        items = [ (10,30), (30,4), (90,9), (13,0), (10, 5), (40, -2), (38, 6), (28,2), (70,9), (11, 3), (25, 4), (99, 100), (12, 21), (32, 23), (78, 87), (65, 56)]
        for key, weight in items:
            skiplist.insert(key, weight)
        if SHAPES and DEBUG: skiplist.print_shape(weightues=True)
        skiplist.remove(10)
        if SHAPES and DEBUG: skiplist.print_shape(weightues=True)
        node = skiplist.find(10)
        self.assertEqual(node.key, -inf)
        skiplist.remove(32)
        node = skiplist.find(32)
        self.assertEqual(node.key, 30)
        skiplist.remove(1000)

