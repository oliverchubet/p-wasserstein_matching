import itertools
import heapq
from math import inf

class PriorityQueue:
    """
    A standard, robust Min-Heap implementation using Python's built-in heapq.
    Implements the exact API expected by ds2 (insert, findmin, remove, changepriority).
    """
    def __init__(self):
        self.pq = []                         # List of entries arranged in a heap
        self.entry_finder = {}               # Mapping of items to entries
        self.counter = itertools.count()     # Unique sequence count to break ties
        self.REMOVED = '<removed-task>'      # Placeholder for removed items

    def insert(self, item, priority):
        """Add a new item or update the priority of an existing item"""
        if item in self.entry_finder:
            self.remove(item)
        count = next(self.counter)
        # Entry structure: [priority, count, item]
        entry = [priority, count, item]
        self.entry_finder[item] = entry
        heapq.heappush(self.pq, entry)

    def remove(self, item):
        """Mark an existing item as removed. Raise KeyError if not found."""
        if item in self.entry_finder:
            entry = self.entry_finder.pop(item)
            entry[-1] = self.REMOVED

    def findmin(self):
        """Return the item with the minimum priority without removing it."""
        # Clean the top of the heap of any removed items
        while self.pq:
            priority, count, item = self.pq[0]
            if item == self.REMOVED:
                heapq.heappop(self.pq)
                continue
            return item
        return None

    def popmin(self):
        """Remove and return the item with the minimum priority."""
        while self.pq:
            priority, count, item = heapq.heappop(self.pq)
            if item != self.REMOVED:
                del self.entry_finder[item]
                return item
        return None

    def changepriority(self, item, priority):
        """Update the priority of an existing item."""
        # The safest way using standard libraries is to remove and re-insert
        self.remove(item)
        self.insert(item, priority)

    def __contains__(self, item):
        return item in self.entry_finder

    def __len__(self):
        return len(self.entry_finder)

    def isEmpty(self):
        return len(self.entry_finder) == 0

    def __str__(self):
            active_items = []
            for priority, count, item in self.pq:
                if item != self.REMOVED:
                    # Change 'item' to 'item.id' here
                    active_items.append((item.id, priority))

            active_items.sort(key=lambda x: x[1])
            return f"PriorityQueue({active_items})"

# -----------------------------------------------------------------------
# MaxLevelHeap (Safe Version)
# -----------------------------------------------------------------------
# This manages a separate Min-Heap for every level, but negates weights 
# so they act like Max-Heaps.
# -----------------------------------------------------------------------
from collections import defaultdict
class MaxLevelHeap:
    def __init__(self):
        self.heaps = defaultdict(PriorityQueue)
        self.levels = set()
        # CACHE: Maps level -> current best item object
        self.max_item_cache = {}

    def _update_cache(self, level):
        """Helper: checking the heap is expensive, so we only do it on updates."""
        if level in self.heaps and not self.heaps[level].isEmpty():
            # This cleans the heap and finds the new best item
            best_item = self.heaps[level].findmin()
            self.max_item_cache[level] = best_item
        else:
            if level in self.max_item_cache:
                del self.max_item_cache[level]

    def insert(self, item, weight, level):
        self.heaps[level].insert(item, -weight)
        self.levels.add(level)
        # Update cache immediately
        self._update_cache(level)

    def remove(self, item):
        for level in list(self.levels):
            if item in self.heaps[level]:
                self.heaps[level].remove(item)

                # Check if the removed item was the current max
                # If so, we MUST scan the heap to find the new max
                if level in self.max_item_cache and self.max_item_cache[level] == item:
                     self._update_cache(level)

                # Clean up empty levels
                if self.heaps[level].isEmpty():
                    del self.heaps[level]
                    self.levels.remove(level)
                    if level in self.max_item_cache:
                        del self.max_item_cache[level]
                return

    def changepriority(self, item, weight):
        for level in self.levels:
            if item in self.heaps[level]:
                self.heaps[level].changepriority(item, -weight)
                # Always update cache on priority change
                self._update_cache(level)
                return

    def findmax(self, level):
        """
        Instant Lookup. No heap logic. No cleaning.
        """
        return self.max_item_cache.get(level)

    def isEmpty(self):
        return len(self.levels) == 0

    # ... (Keep your __str__ and __contains__ as they were)

#class MaxLevelHeap:
#    def __init__(self):
#        self.heaps = defaultdict(PriorityQueue)
#        self.levels = set()
#
#    def insert(self, item, weight, level):
#        # Negate weight so MinHeap behaves like MaxHeap
#        self.heaps[level].insert(item, -weight)
#        self.levels.add(level)
#
#    def remove(self, item):
#        for level in list(self.levels):
#            if item in self.heaps[level]:
#                self.heaps[level].remove(item)
#                if self.heaps[level].isEmpty():
#                    del self.heaps[level]
#                    self.levels.remove(level)
#                return
#
#    def changepriority(self, item, weight):
#        for level in self.levels:
#            if item in self.heaps[level]:
#                self.heaps[level].changepriority(item, -weight)
#                return
#
#    def findmax(self, level):
#        if level in self.heaps and not self.heaps[level].isEmpty():
#            return self.heaps[level].findmin()
#        return None
#
#    def isEmpty(self):
#        """Returns True if there are no items in any level."""
#        return len(self.levels) == 0

    def __contains__(self, item):
        for h in self.heaps.values():
            if item in h: return True
        return False

    def __str__(self):
            if not self.levels:
                return "MaxLevelHeap(Empty)"

            lines = ["MaxLevelHeap:"]
            for level in sorted(self.levels):
                pq = self.heaps[level]
                active_items = []

                for priority, count, item in pq.pq:
                    if item != pq.REMOVED:
                        # Change 'item' to 'item.id' here
                        active_items.append((item.id, -priority))

                active_items.sort(key=lambda x: x[1], reverse=True)
                lines.append(f"  Level {level}: {active_items}")

            return "\n".join(lines)
