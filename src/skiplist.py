from math import inf
from random import choice

TRACE = True

class SkipList(dict):
    def __init__(self, height=2):
        dict.__init__(self)
        self.head = Node(-inf, -inf)
        self.head.next = Node(inf, -inf)
        self.height = 1
        while self.height < height:
            self._increase_height()

    def insert(self, key, weight):
        if TRACE: print("SkipList.insert(", key, ",", weight, ")")
        h = 1
        while choice([True, False]) and h < self.height: h+=1
        if h == self.height: self._increase_height()
        parent = self.head
        parent = self.find(key, depth=self.height-h)
        node = Node(key, weight)
        self[key] = node
        parent.connect(node)
        while parent.down is not None:
            parent = self.find(key, depth=1, handle=parent)
            node.down = Node(key, weight)
            node = node.down
            parent.connect(node)

    def remove(self, key):
        if TRACE: print("SkipList.remove(", key, ")")
        if key not in self: return
        cur = self.head
        while cur.down is not None:
            cur = self.find(key, depth=1, handle=cur)
            while cur.next.key == key:
                cur.next = cur.next.next
        return self.pop(key, None)

    def find(self, key, depth=None, handle=None):
        if TRACE: print("SkipList.find(", key, "depth=", depth, "handle=", handle, ")")
        cur = handle
        if cur is None: cur = self.head
        if depth is None: depth = self.height
        while True:
            while cur.next.key < key:
                cur = cur.next
            if depth > 0 and cur.down is not None:
                cur = cur.down
                depth -= 1
            else: return cur

    def print_shape(self, keys=False, weightues=False):
        if TRACE: print("SkipList.print_shape(keys=", keys, "weightues=", weightues, ")")
        #print("SkipList.print_shape(keys=", keys, "weightues=", weightues, ")")
        print("SkipList Shape:")
        handle = self.head
        while handle is not None:
            left, right = handle, handle.next
            while right is not None:
                if weightues: 
                    if left.weight == -inf:
                        print("-o", end=" ")
                    else:
                        print(left.weight, end=" ")
                elif keys: print(left.key, end=" ")
                else: print("o-", end="")
                lcur, rcur = left, right
                while lcur.down is not None:
                    lcur = lcur.down
                while rcur.down is not None:
                    rcur = rcur.down
                while lcur.next is not rcur:
                    lcur = lcur.next
                    if weightues or keys: print("--", end=" ")
                    else: print("--", end="")
                left, right = right, right.next
            print("o")
            if handle.down is None:
                left, right = handle, handle.next
                while right is not None:
                    print("--", end=" ")
                    left, right = right, right.next
                print()
                left, right = handle, handle.next
                while right is not None:
                    if left.key == -inf: print("-o", end=" ")
                    else: print(left.key, end=" ")
                    left, right = right, right.next
                print("\n\n")
            handle = handle.down

    def _increase_height(self):
        self.height += 1
        new_head = Node(self.head.key, self.head.weight)
        new_head.down = self.head
        new_tail = Node(inf, -inf)
        new_tail.down = self.head.next
        new_head.next = new_tail
        self.head = new_head


class Node:
    def __init__(self, key, weight):
        self.key = key
        self.max_key = key
        self.weight = weight 
        self.next = None
        self.down = None

    def connect(self, other):
        other.next = self.next
        self.next = other

    def __str__(self):
        #return "".join(["(", str(self.key), ",", str(self.weight),")"]) 
        return "".join(["(", str(self.key), ",", str(self.weight),",", str(self.max_key),")"]) 
