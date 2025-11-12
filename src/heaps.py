from ds2.priorityqueue import PriorityQueue

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

    def __contains__(self, item):
        return item in self._itemmap

class MinHeap(PriorityQueue):

    def __contains__(self, item):
        return item in self._itemmap

