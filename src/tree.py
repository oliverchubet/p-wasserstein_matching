class Tree:
    def __init__(self):
        self.root = Node("root")
        self.nodes = {"root": self.root}

    def connect(self, point, parent):
        parent_node = self.nodes[parent]
        node = Node(point, parent_node)
        self.nodes[point] = node
        parent_node.children.add(point)

    def remove(self, point):
        node = self.nodes[point]
        parent = node.parent
        for child in node.children:
            self.remove(child)
        if parent is not None:
            return parent.pt

    def path(self, point):
        cur = self.nodes[point]
        path = [point]
        while cur is not self.root:
            cur = cur.parent
            path.append(cur.pt)
        return path


class Node:
    def __init__(self, pt, parent=None):
        self.pt = pt
        self.parent = parent
        self.children = set()

