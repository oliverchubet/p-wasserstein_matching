import unittest
from skipheap import SkipHeap
from math import inf
from random import randint

SHAPES = False
DEBUG = False

class TestSkipHeap(unittest.TestCase):

    def test_skipheap_init(self):
        skeap = SkipHeap()
        node = skeap.find(-inf)
        self.assertEqual(node.weight, -inf)
        node = skeap.find(0)
        self.assertEqual(node.weight, -inf)
        node = skeap.find(1000)
        self.assertEqual(node.weight, -inf)

    def test_skipheap_increase_height(self):
        skeap = SkipHeap()
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
        if SHAPES: skeap.print_shape()

    def test_skipheap_insert(self):
        skeap = SkipHeap()
        n = 100
        if SHAPES and DEBUG:
            n = 25
            skeap.print_shape()
        for i in range(n):
            key, weight = randint(10,99), randint(10,99)
            #print(i, ";", key, weight, end=" ")
            if key not in skeap:
                skeap.insert(key, weight)
            start = skeap.head
            cur = start
            #for j in range(skeap.height):
            #    while cur.next is not None:
            #        print(cur.key, end=" ")
            #        cur = cur.next
            #    print()
            #    start = start.down
            #    cur = start
            if SHAPES and DEBUG: skeap.print_shape(weightues=True)
        if SHAPES and not DEBUG: skeap.print_shape()

    def test_skipheap_incremental_insert(self):
        skeap = SkipHeap()
        if SHAPES: skeap.print_shape()
        for i in range(10, 35):
            key, weight = i, randint(10,99)
            if key not in skeap:
                skeap.insert(key, weight)
            if SHAPES: skeap.print_shape()

    def test_skipheap_keys_in_order(self):
        skeap = SkipHeap()
        for i in range(100):
            key, weight = i, randint(10,99)
            if key not in skeap:
                skeap.insert(key, weight)
        start = skeap.head
        while start is not None:
            left, right = start, start.next
            while right is not None:
                self.assertTrue(left.key < right.key)
                left, right = right, right.next
            start = start.down

    def test_skipheap_weights(self):
        skeap = SkipHeap()
        for i in range(25):
            key, weight = i, randint(10,99)
            if key not in skeap:
                skeap.insert(key, weight)
                if SHAPES: skeap.print_shape(weightues=True)
        start = skeap.head
        while start is not None:
            left, right = start, start.next
            while right is not None:
                lcur, rcur = left, right
                while lcur.down is not None:
                    lcur = lcur.down
                while rcur.down is not None:
                    rcur = rcur.down
                while lcur.next is not rcur:
                    lcur = lcur.next
                    if lcur.weight > left.weight:
                        print("lcur.weight > left.weight")
                        skeap.print_shape(weightues=True)
                    self.assertTrue(lcur.weight <= left.weight)
                left, right = right, right.next
            start = start.down
        for i in range(10):
            skeap.remove(i)
            #if SHAPES: skeap.print_shape(weightues=True)
            skeap.print_shape(weightues=True)
        start = skeap.head
        while start is not None:
            left, right = start, start.next
            while right is not None:
                lcur, rcur = left, right
                while lcur.down is not None:
                    lcur = lcur.down
                while rcur.down is not None:
                    rcur = rcur.down
                while lcur.next is not rcur:
                    lcur = lcur.next
                    if lcur.weight > left.weight:
                        print("lcur.key =", lcur.key, "; left.key =", left.key)
                        print("lcur.weight =",lcur.weight, "> left.weight =", left.weight)
                        skeap.print_shape(weightues=True)
                    self.assertTrue(lcur.weight <= left.weight)
                left, right = right, right.next
            start = start.down

    def test_skipheap_update(self):
        skeap = SkipHeap()
        skeap.insert(1,2)
        skeap.insert(8,1)
        skeap.insert(35,9)
        node = skeap.findmax(30)
        self.assertEqual(node.weight, 2)
        node = skeap.findmax(36)
        self.assertEqual(node.weight, 9)
        skeap.update(35, 10)
        node = skeap.findmax(36)
        self.assertEqual(node.weight, 10)
        node = skeap.findmax(30)
        self.assertEqual(node.weight, 2)
        start = skeap.head
        while start is not None:
            left, right = start, start.next
            while right is not None:
                lcur, rcur = left, right
                while lcur.down is not None:
                    lcur = lcur.down
                while rcur.down is not None:
                    rcur = rcur.down
                while lcur.next is not rcur:
                    lcur = lcur.next
                    if lcur.weight > left.weight:
                        print("lcur.weight > left.weight")
                        skeap.print_shape(weightues=True)
                    self.assertTrue(lcur.weight <= left.weight)
                left, right = right, right.next
            start = start.down


    def test_skipheap_remove(self):
        skeap = SkipHeap()
        items = [ (30,14), (90,-9), (13,-1), (10, 15), (40, -2), (38, 26), (28,32), (70,19), (11, -3), (25, 14), (99, 91), (12, 21), (32, 23), (78, 87), (65, 56)]
        for key, weight in items:
            skeap.insert(key, weight)
        #if SHAPES and DEBUG: skiplist.print_shape(weightues=True)
        if SHAPES: skeap.print_shape(weightues=True)
        skeap.remove(10)
        if SHAPES: skeap.print_shape(weightues=True)
        #if SHAPES and DEBUG: skeap.print_shape(weightues=True)
        node = skeap.find(10)
        self.assertEqual(node.key, -inf)
        skeap.remove(32)
        if SHAPES: skeap.print_shape(weightues=True)
        node = skeap.find(32)
        self.assertEqual(node.key, 30)
        skeap.remove(1000)


    def test_skipheap_findmax(self):
        skeap = SkipHeap()

        skeap.insert(1, 2)

        #skeap.print_shape(weightues=True)

        self.assertEqual(len(skeap), 1)
        node = skeap.findmax(3)
        self.assertEqual(node.weight, 2)
        node = skeap.findmax(1)
        self.assertEqual(node.weight, -inf)
        node = skeap.findmax(0)
        self.assertEqual(node.weight, -inf)

        skeap.insert(5,10)

        #skeap.print_shape(weightues=True)

        self.assertEqual(len(skeap), 2)
        node = skeap.findmax(1)
        self.assertEqual(node.weight, -inf)
        node = skeap.findmax(2)
        self.assertEqual(node.weight, 2)
        node = skeap.findmax(3)
        self.assertEqual(node.weight, 2)
        node = skeap.findmax(5)
        self.assertEqual(node.weight, 2)
        node = skeap.findmax(6)
        self.assertEqual(node.weight, 10)
        node = skeap.findmax(1000)
        self.assertEqual(node.weight, 10)
