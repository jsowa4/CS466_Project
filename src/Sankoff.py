import Node
import sys

def cost(s,t):
    """
    Computes the cost of a mutation between s and t
    """
    return 0 if s==t else 1


def sankoff_fill(root_node, states_alphabet):
    """
    Takes a reference to the root of a character based phylogeny tree,
    and fills in the mu values for all nodes

    """
    if root_node.leaf:
        for s in states_alphabet:
            if root_node.state == s:
                root_node.mu[s] = 0
            else:
                root_node.mu[s] = float("inf")
    else:
        # recurse (post order traversal)
        for w in root_node.children:
            sankoff_fill(w, states_alphabet)

        # post order processing of current node
        for s in states_alphabet:
            root_node.mu[s] = 0.0
            for w in root_node.children:
                min_val = float("inf")
                for t in states_alphabet:
                    if min_val > cost(s,t) + w.mu[t]:
                        min_val = cost(s,t) + w.mu[t]
                root_node.mu[s] += min_val

def sankoff_backtrace(root_node, state_alphabet):
    """
    Takes a reference to the root node of a character based phylogeny tree and
    fills in the state values with the appropriate character.

    """
    # base case
    if root_node.leaf:
        return

    # process current node
    if root_node.root:
        min_val = float("inf")
        argmin_val = None
        for s in root_node.mu:
            if root_node.mu[s] < min_val:
                min_val = root_node.mu[s]
                argmin_val = s
        root_node.state = argmin_val
        root_node.min_val = min_val
    else:
        s = root_node.parent.state
        min_val = float("inf")
        min_t = None
        for t in state_alphabet:
            if cost(s,t) + root_node.mu[t] < min_val:
                min_val = cost(s,t) + root_node.mu[t]
                min_t = t
        root_node.state = min_t
        root_node.min_val = min_val

    # recurse
    for w in root_node.children:
        sankoff_backtrace(w, state_alphabet)


# Parse the tree structure from the file
def build_tree(tree_file, labels_file, states_alphabet):
    # Load node labels and states
    with open(labels_file, 'r') as f:
        labels_data = f.readlines()
    node_labels = {}
    for line in labels_data:
        parts = line.strip().split('\t')
        node_id = parts[0]
        node_states = parts[1:]  # States are in the rest of the line
        node_labels[node_id] = node_states

    # Create a dictionary of nodes
    nodes = {}
    for node_id, states in node_labels.items():
        is_leaf = len(states) > 0
        # Assign the first state as the current character; this can be updated later
        character = states[0] if is_leaf else None
        nodes[node_id] = Node(node_id, character, states_alphabet, leaf=is_leaf)

    # Build the tree structure from the tree file
    with open(tree_file, 'r') as f:
        tree_data = f.readlines()
    for line in tree_data:
        parent_id, child_id = line.strip().split()
        parent_node = nodes[parent_id]
        child_node = nodes[child_id]
        parent_node.add_child(child_node)

    # Identify the root (a node without a parent)
    root = None
    for node in nodes.values():
        if node.parent is None:
            node.root = True
            root = node
            break

    return root, nodes

def main():
    states_alphabet = ['LdfBx', 'LdfRx', 'RIIM0', 'RulM0']
    tree_file = '../data/osteosarcoma/OSCE2.tree'
    labels_file = '../data/osteosarcoma/OSCE2.observed.labeling'

    root, nodes = build_tree(tree_file, labels_file, states_alphabet)

    # Verify the structure
    print("Root:", root.name)
    for node_name, node in nodes.items():
        print(f"Node {node_name}: Parent={node.parent.name if node.parent else 'None'}, "
            f"Children={[child.name for child in node.children]}")