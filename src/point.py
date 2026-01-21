
class Point:

    count = 0

    def __init__(self, coords, weight):
        self.id = Point.count
        Point.count += 1
        self.coords = coords 
        self.dual_weight = 0
        self.len_path = None
        self.match = None
        self.inserted = False

    #def __eq__(self, other):
    #    return self.id == other.id

    def __str__(self):
        return "p" + str(self.id) #+": "+str(self.coords)

class Edge:
    def __init__(self, a, b, slack, center, level):
        self.a = a
        self.b = b
        self.slack = slack
        self.cluster = center
        self.level = level

    def __str__(self):
        if self.a is None: str_a = "None"
        else: str_a = str(self.a.id)
        if self.b is None: str_b = "None"
        else: str_b = str(self.b.id)
        return "".join(["Edge (", str_a, ",", str_b, "); slack = ", str(self.slack),"; ctr =", str(self.cluster.id)])
