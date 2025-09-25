
class Point:

    count = 0

    def __init__(self, coords, weight):
        self.id = Point.count
        Point.count += 1
        self.coords = coords 
        self.dual_weight = 0
        self.weight = None
        self.len_path = None
        self.free = True
        self.match = None
        self.matched = False
        self.visited = False
        self.DFS_visited = False

class Edge:
    def __init__(self, a, b, slack, center, level):
        self.a = a
        self.b = b
        self.slack = slack
        self.cluster = center # the cluster that created it
        self.level = level
