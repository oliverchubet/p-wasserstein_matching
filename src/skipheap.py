from math import inf
from random import choice
from skiplist import SkipList

TRACE = True

class SkipHeap(SkipList):

    def insert(self, key, weight):
        if TRACE: print("SkipHeap.insert(", key, ",", weight, ")")
        assert(key not in self)
        SkipList.insert(self, key, weight)

        handle = self.head
        while handle.down is not None:
            if handle.next.key > key:
                if handle.weight < weight:
                    handle.weight, handle.max_key = weight, key
            else:
                if handle.max_key > key:
                    if handle.weight > weight:
                        handle.next.weight, handle.next.max_key = handle.weight, handle.max_key
                    node = self.findmax(key, handle=handle)
                    handle.weight, handle.max_key = node.weight, node.max_key
            handle = self.find(key, depth=1, handle=handle)

        handle = self[key]
        while handle.down is not None:
            if handle.max_key == key:
                node = self._openfindmax(key, handle.next.key, handle)
                if node.weight > weight:
                    handle.weight, handle.max_key = node.weight, node.max_key
                else: break
            handle = handle.down
                


    def _increase_height(self):
        if TRACE: print("SkipHeap._increase_height()")
        SkipList._increase_height(self)
        self.head.max_key = self.head.down.max_key

    def update(self, key, weight):
        if TRACE: print("SkipHeap.update(", key, ",", weight, ")")
        handle = self.head
        while handle.down is not None:
            if handle.weight < weight:
                handle.weight = weight
                handle.max_key = key
            elif handle.weight > weight and handle.max_key == key:
                left = self.findmax(key, handle=handle)
                print("left", left)
                #right = self.openfindmax(self[key], end=handle.next)
                right = self._openfindmax(key, handle.next.key, handle=self[key])
                print("right", right)
                if left is None: #and left.weight >= weight:
                    handle.weight = weight
                    handle.max_key = key
                else:
                    handle.weight = left.weight
                    handle.max_key = left.max_key
                if right is not None and right.weight >= handle.weight:
                    handle.weight = right.weight
                    handle.max_key = right.max_key
            #if handle.key < key:
            if handle.down is None: break
            handle = self.find(key, depth=1, handle=handle)
            #elif handle.next.key == key:
            if handle.next.key == key:
                handle = handle.next
            #else: handle = handle.down
        handle.weight = weight

    def remove(self, key):
        if TRACE: print("SkipHeap.remove(", key, ")")
        self.update(key, -inf)
        self.print_shape(weightues=True)
        SkipList.remove(self, key)
        #if key not in self: return
        #SkipList.remove(self, key)
        #handle = self.head
        #while handle is not None:
        #    node = self.findmax(key, handle=handle)
        #    handle.weight, handle.max_key = node.weight, node.max_key
        #    node = self._openfindmax(key, handle.next.key, handle=handle)
        #    if node.weight >= handle.weight:
        #        handle.weight, handle.max_key = node.weight, node.max_key
        #    if handle.down is None: break
        #    handle = self.find(key, depth=1)

    # max over the half open interweight [handle.key, key)
    def findmax(self, key, handle=None):
        if TRACE: print("SkipHeap.findmax(", key, ", handle=", handle, ")") 
        if handle is None: handle = self.head
        max_found, max_weight = None, -inf
        while handle is not None:
            while handle.next.key < key:
                if handle.weight > max_weight:
                    max_found, max_weight = handle, handle.weight
                handle = handle.next
            if handle.max_key < key:
                if handle.weight >= max_weight:
                    max_found, max_weight = handle, handle.weight
                return max_found
            handle = handle.down
        return max_found
            #print("handle", handle)
            #if handle.max_key < key and max_weight <= handle.weight:
            #    max_found, max_weight = handle, handle.weight
            #    if handle.next.key >= key:
            #        return max_found
            #if handle.next.key < key: handle = handle.next
            #else: handle = handle.down
        #return max_found

    # max over the open interweight (handle.key, key)
    #def _openfindmax(self, handle, key):
    def _openfindmax(self, left_key, right_key, handle):
        #if TRACE: print("SkipHeap.openfindmax(handle=", handle, ", end=", end, ")") 
        if TRACE: print("SkipHeap._openfindmax(handle=", handle, ", left_key =", left_key, ", right_key =", right_key, ")") 
        while handle.max_key < left_key:
            if handle.next.key > left_key:
                if handle.down is None: break
                handle = handle.down
            else:
                handle = handle.next
        if handle.down is None and handle.key == left_key:
            if handle.key == right_key: return None
            return self.findmax(right_key, self[handle.next.key])
        else:
            max_found, max_weight = handle, handle.weight
            if handle.key == -inf:
                found = self.findmax(right_key)
            else:
                found = self.findmax(right_key, self[handle.key])
            if found.weight > max_weight: return found
            else: return max_found
        #max_found, max_weight = None, -inf
        #while handle.next.key > end.key:
        #    handle = handle.down
        #while handle.down is not None:
        #    if handle.max_key > handle.key:
        #        if max_weight <= handle.weight:
        #            return self[handle.max_key]
        #    handle = handle.down
        #    node = self.findmax(end.key, handle=handle.next)
        #    if node is not None and node.weight > max_weight:
        #        max_found, max_weight = node, node.weight

