class Node:
    def __init__(self, name,  character, states_alphabet, root=False, leaf=False):
        self.name = name
        self.state = character
        self.leaf = leaf
        self.root = root
        self.parent = None
        self.children = []
        self.mu = {s:float("inf") for s in states_alphabet}
        self.min_val = None

    def add_child(self, child):
        self.children.append(child)
        child.assign_parent(self)

    def assign_parent(self, parent):
        self.parent = parent