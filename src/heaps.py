from ds2.priorityqueue import PriorityQueue
from collections import defaultdict

class MaxLevelHeap:
    def __init__(self):
        self.heaps = defaultdict(MaxHeap)
        self._itemmap = dict()
        self.levels = set()
        self.max_item = None
        self.max_val = None

    def findmax(self, level=1):
        self.max_item = None
        self.max_val = None
        for i in self.levels:
            if i is None or (level is not None and i <= level):
                item = self.heaps[i].findmax()
                val = -self.heaps[i].findmaxval()
                if self.max_val is None or val > self.max_val:
                    self.max_item = item
                    self.max_val = val
                    #if self.max_val is not None:
                        #if val == -10: print(i, self.max_item.id)
        #if self.max_item is not None:
        #    print(self.max_item.id, self.max_val)
        #    print("++++++++++++")
        return self.max_item

    def findmaxval(self, level):
        if len(self.heaps[level]) > 0:
            return self.heaps[level].findmaxval()

    def insert(self, item, priority, level):
        self.heaps[level].insert(item, priority)
        self._itemmap[item] = level
        self.levels.add(level)

    def remove(self, item):
        level = self._itemmap.pop(item, None)
        if item in self.heaps[level]:
            self.heaps[level].remove(item)
        if len(self.heaps[level]) == 0 and level in self.levels:
            self.levels.remove(level)

    def changepriority(self, item, priority):
        level = self._itemmap[item]
        self.heaps[level].changepriority(item, priority)

    def isEmpty(self):
        return len(self.levels) == 0

    def __contains__(self, item):
        return item in self._itemmap

    def __str__(self):
        self.findmax()
        s = ""
        s += "max_val " + str(self.max_val) + " "
        if self.max_item is not None:
            s += "max_item " + str(self.max_item.id)
        s += "\n"
        for level in self.levels:
            str_len = 0
            s += "[" + str(level)+ ":"
            for entry in self.heaps[level]._entries:
                s += "".join(["(",str(entry.item.id),",", str(entry.priority),")"]) 
            s += "]\n"
        return s

class MinLevelHeap(MaxLevelHeap):
    def insert(self, item, priority, level):
        MaxLevelHeap.insert(self, item, -priority, level)

    def changepriority(self, item, priority):
        MaxLevelHeap.insert(self, item, -priority)


class MaxHeap(PriorityQueue):
    def __init__(self):
        super().__init__()
        self.removemax = self.removemin
        self.findmax = self.findmin
        self.len = 0

    def insert(self, item, priority):
        PriorityQueue.insert(self, item, -priority)
        self.len += 1

    def changepriority(self, item, priority):
        PriorityQueue.changepriority(self, item, -priority)

    def findmaxval(self):
        if self.len > 0:
            return self._entries[0].priority

    #def findmax(self):
    #    if len(self) > 0:
    #        return self._entries[0].item

    def __contains__(self, item):
        return item in self._itemmap

class MinHeap(PriorityQueue):

    def __contains__(self, item):
        return item in self._itemmap

    def __str__(self):
        s = ""
        for entry in self._entries:
            s += "".join(["(",str(entry.item.center.id),",",str(entry.item.min_slack_edge),", priority =", str(entry.priority),")"]) + "\n"
        return s
