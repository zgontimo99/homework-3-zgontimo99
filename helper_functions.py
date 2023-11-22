"""Helper functions for HW3"""
import numpy as np
from copy import deepcopy
from matplotlib.axes import Axes


class Node:
    def __init__(
        self,
        name: str,
        left: "Node",
        left_distance: float,
        right: "Node",
        right_distance: float,
        confidence: float = None,
        count: int = 0,
        lefts: int = 0,
        y_pos: float = 0,
    ):
        """A node in a binary tree produced by neighbor joining algorithm.

        Parameters
        ----------
        name: str
            Name of the node.
        left: Node
            Left child.
        left_distance: float
            The distance to the left child.
        right: Node
            Right child.
        right_distance: float
            The distance to the right child.
        confidence: float
            The confidence level of the split determined by the bootstrap method.
            Only used if you implement Bonus Problem 1.

        Notes
        -----
        The current public API needs to remain as it is, i.e., don't change the
        names of the properties in the template, as the tests expect this kind
        of structure. However, feel free to add any methods/properties/attributes
        that you might need in your tree construction.

        """
        self.name = name
        self.left = left
        self.left_distance = left_distance
        self.right = right
        self.right_distance = right_distance
        self.confidence = confidence
        self.count = count
        self.lefts = lefts
        self.y_pos = y_pos



def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("the brown cat", "these brownies", lambda x, y: [-1, 1][x == y])
    ('the-- brown cat', 'these brownies-', 3.0)

    Other alignments are also possible.

    """
    n = len(seq1)
    m = len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + scoring_function(seq1[i - 1],'*')
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + scoring_function('*',seq2[j - 1])

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                match_score = dp[i - 1][j - 1] + scoring_function(seq1[i - 1],seq2[j - 1])
            else:
                match_score = dp[i - 1][j - 1] + scoring_function(seq1[i - 1],seq2[j - 1])
            gap1_score = dp[i - 1][j] + scoring_function(seq1[i - 1],'*')
            gap2_score = dp[i][j - 1] + scoring_function('*',seq2[j - 1])
            dp[i][j] = max(match_score, gap1_score, gap2_score)

    
    align1, align2 = '', ''
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and dp[i][j] == dp[i - 1][j] + scoring_function(seq1[i - 1],'*'):
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + scoring_function('*',seq2[j - 1]):
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1
        else:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1

    return align1, align2, dp[n][m]





def neighbor_joining(distances: np.ndarray, labels: list) -> Node:
    """The Neighbor-Joining algorithm.

    For the same results as in the later test dendrograms;
    add new nodes to the end of the list/matrix and
    in case of ties, use np.argmin to choose the joining pair.

    Parameters
    ----------
    distances: np.ndarray
        A 2d square, symmetric distance matrix containing distances between
        data points. The diagonal entries should always be zero; d(x, x) = 0.
    labels: list
        A list of labels corresponding to entries in the distances matrix.
        Use them to set names of nodes.

    Returns
    -------
    Node
        A root node of the neighbor joining tree.

    """
    n = len(distances)
    nodes={}
    m=0
    for label in labels:
        nodes[label]=Node(label, None, 0, None, 0, y_pos=n-m)
        m=m+1
    
    while n > 2:
        Q = np.zeros_like(distances)
        total_distance = np.sum(distances, axis=1)

        for i in range(n):
            for j in range(i + 1, n):
                Q[i, j] = (n - 2) * distances[i, j] - total_distance[i] - total_distance[j]
                Q[j, i] = Q[i, j]
        agh=np.argmin(Q)
        min_i=int(agh/n)
        min_j=agh%n
        print(Q)
        limb_i = 0.5 * (distances[min_i, min_j] + (total_distance[min_i] - total_distance[min_j]) / (n - 2))
        limb_j = distances[min_i, min_j] - limb_i
        print(limb_i, limb_j)
        X=Node(labels[min_i]+labels[min_j],nodes[labels[min_i]],limb_i,nodes[labels[min_j]],limb_j)
        nodes[labels[min_i]+labels[min_j]]=X
        nodes.pop(labels[min_i])
        nodes.pop(labels[min_j])
        limb_i = 0.5 * (distances[min_i, min_j] + (total_distance[min_i] - total_distance[min_j]) / (n - 2))
        labels_new=[]
        for i in range(n):
            if i == min_i:
                labels_new.append(labels[min_i]+labels[min_j])
            elif i != min_j:
                labels_new.append(labels[i])
        
        new_distances = np.zeros((n - 1, n - 1))

        for i in range(n):
            m=i-1
            if i != min_i and i != min_j and min_j>i:
                new_distances[min_i, i] = (distances[min_i, i] + distances[min_j, i] - distances[min_i,min_j]) / 2

                new_distances[i, min_i] = new_distances[min_i, i]
            elif i != min_i and i != min_j:
                new_distances[min_i, m] = (distances[min_i, i] + distances[min_j, i] - distances[min_i,min_j]) / 2

                new_distances[m, min_i] = new_distances[min_i, m]

        for i in range(n):
            for j in range(n):
                if i != min_i and j != min_i and i != min_j and j != min_j and min_j>i:
                    new_distances[i, j] = distances[i, j]
                elif i != min_i and j != min_i and i != min_j and j != min_j:
                    new_distances[i-1, j-1] = distances[i, j]
        print(new_distances)
        distances=new_distances
        labels=labels_new
        print(labels)
        n -= 1

    X=Node(labels[0]+labels[1],nodes[labels[0]],distances[0,1]/2,nodes[labels[1]],distances[0,1]/2)
    nodes[labels[0]+labels[1]]=X
    nodes.pop(labels[0])
    nodes.pop(labels[1])
    for element in nodes.values():
        Y=element
    
    return Y
    

def plot_nj_tree(tree: Node, ax: Axes = None, **kwargs) -> None:
    """A function for plotting neighbor joining phylogeny dendrogram.

    Parameters
    ----------
    tree: Node
        The root of the phylogenetic tree produced by `neighbor_joining(...)`.
    ax: Axes
        A matplotlib Axes object which should be used for plotting.
    kwargs
        Feel free to replace/use these with any additional arguments you need.
        But make sure your function can work without them, for testing purposes.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>>
    >>> tree = neighbor_joining(distances)
    >>> fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    >>> plot_nj_tree(tree=tree, ax=ax)
    >>> fig.savefig("example.png")

    """
    def get_leaf_count(tree: Node):
        tree.count=0
        if tree.left!=None:
            tree.count+=get_leaf_count(tree.left)
        if tree.right!=None:
            tree.count+=get_leaf_count(tree.right)
        if tree.count == 0:
            tree.count = 1
        return tree.count
    def get_left_count(tree: Node):
        if tree is None:
            return 0
        left_leaves = get_left_count(tree.left)
        right_leaves = get_left_count(tree.right)
        if tree.left is None and tree.right is None:
            left_leaves += 1
        tree.lefts = left_leaves
        return left_leaves + right_leaves
        return left_leaves + right_leaves, left_count + right_count
        
    def get_node_ys(tree: Node,n):
        tree.y_pos=tree.left.count/tree.count*(tree.count-1)+tree.left.lefts
        if tree.left.count > 1 and tree.left.left != None:
            get_node_ys(tree.left,n)
        if tree.right.count > 1 and tree.right.left != None:
            get_node_ys(tree.right,n)
    def plot_node(tree:Node, val):
        ax.plot((val,val+tree.left_distance),(tree.left.y_pos,tree.left.y_pos),color='black')
        ax.plot((val,val+tree.right_distance),(tree.right.y_pos,tree.right.y_pos),color='black')
        ax.vlines(val,tree.left.y_pos,tree.right.y_pos,color='black')
        if tree.left.count>1 and tree.left.left != None:
            plot_node(tree.left,val+tree.left_distance)
        if tree.right.count>1 and tree.right.left != None:
            plot_node(tree.right,val+tree.right_distance)

    get_leaf_count(tree)
    get_left_count(tree)
    n=tree.count
    get_node_ys(tree,n)
    plot_node(tree, 0)
    plt.show()
    return ax


def _find_a_parent_to_node(tree: Node, node: Node) -> tuple:
    """Utility function for reroot_tree"""
    stack = [tree]

    while len(stack) > 0:

        current_node = stack.pop()
        if node.name == current_node.left.name:
            return current_node, "left"
        elif node.name == current_node.right.name:
            return current_node, "right"

        stack += [
            n for n in [current_node.left, current_node.right] if n.left is not None
        ]

    return None


def _remove_child_from_parent(parent_node: Node, child_location: str) -> None:
    """Utility function for reroot_tree"""
    setattr(parent_node, child_location, None)
    setattr(parent_node, f"{child_location}_distance", 0.0)


def reroot_tree(original_tree: Node, outgroup_node: Node) -> Node:
    """A function to create a new root and invert a tree accordingly.

    This function reroots tree with nodes in original format. If you
    added any other relational parameters to your nodes, these parameters
    will not be inverted! You can modify this implementation or create
    additional functions to fix them.

    Parameters
    ----------
    original_tree: Node
        A root node of the original tree.
    outgroup_node: Node
        A Node to set as an outgroup (already included in a tree).
        Find it by it's name and then use it as parameter.

    Returns
    -------
    Node
        Inverted tree with a new root node.
    """
    tree = deepcopy(original_tree)

    parent, child_loc = _find_a_parent_to_node(tree, outgroup_node)
    distance = getattr(parent, f"{child_loc}_distance")
    _remove_child_from_parent(parent, child_loc)

    new_root = Node("new_root", parent, distance / 2, outgroup_node, distance / 2)
    child = parent

    while tree != child:
        parent, child_loc = _find_a_parent_to_node(tree, child)

        distance = getattr(parent, f"{child_loc}_distance")
        _remove_child_from_parent(parent, child_loc)

        empty_side = "left" if child.left is None else "right"
        setattr(child, f"{empty_side}_distance", distance)
        setattr(child, empty_side, parent)

        if tree.name == parent.name:
            break
        child = parent

    other_child_loc = "right" if child_loc == "left" else "left"
    other_child_distance = getattr(parent, f"{other_child_loc}_distance")

    setattr(child, f"{empty_side}_distance", other_child_distance + distance)
    setattr(child, empty_side, getattr(parent, other_child_loc))

    return new_root


def sort_children(tree: Node) -> None:
    """Sort the children of a tree by their corresponding number of leaves.

    The tree can be changed inplace.

    Paramteres
    ----------
    tree: Node
        The root node of the tree.

    """
    raise NotImplementedError()


def plot_nj_tree_radial(tree: Node, ax: Axes = None, **kwargs) -> None:
    """A function for plotting neighbor joining phylogeny dendrogram
    with a radial layout.

    Parameters
    ----------
    tree: Node
        The root of the phylogenetic tree produced by `neighbor_joining(...)`.
    ax: Axes
        A matplotlib Axes object which should be used for plotting.
    kwargs
        Feel free to replace/use these with any additional arguments you need.

    Example
    -------
    >>> import matplotlib.pyplot as plt
    >>>
    >>> tree = neighbor_joining(distances)
    >>> fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))
    >>> plot_nj_tree_radial(tree=tree, ax=ax)
    >>> fig.savefig("example_radial.png")

    """
    raise NotImplementedError()
