"""
Attempt to code network propagation
"""


import networkx as nx
import numpy as np
import os
# from stringInteractions2namedInteractions import create_aliasdict
#
# seedlist = r"txt/string_seeds.txt"
# aliasfile = "gz/39947.protein.aliases.v11.5.txt.gz"
regen = False

# output numpy file name (without the .npy)
outfile = "p"

scale_limit = 100.0


def weights_from_seeds(graph, seedlist, weight=100):
    """
    Takes an input graph and creates a list of weights for the graph nodes, where each seed node from a list
    of seeds will get assigned the starting weight value (defaults to 100) and all other nodes will start at
    0

    :param graph: input graph
    :param seedlist: list of seed nodes to assign initial weight to
    :param weight: starting weight value
    :return: List of weights corresponding to each node in the graph
    """
    p0 = []
    # Get seed proteins and assign to them an input weight of 100 and 0 to the rest
    # with open(os.path.join(seedlist), 'r') as file:
    #     seeds = [line.rstrip() for line in file.readlines()]

    seeds= seedlist

    for name in graph.nodes:
        if name in seeds:
            p0.append(weight)
        else:
            p0.append(0)

    p0 = np.array(p0)


    return p0


def _rwr(p0, alpha, A, invD, iter):
    """
    Random walk with restart implementation that takes in initialization parameters and only carries out
    the walk. The initialization parameters have been taken out of the method to allow relatively constant elements
    of the calculation to be provided once without the need to recalculate multiple times. Implements cupy to
    speed up the walk steps

    :param p0: initial weights
    :param alpha: learning rate
    :param A: adjacency matrix of graph
    :param invD: inverse of the degree matrix of the graph
    :param iter: number of iterations to restart the walk
    :return: numpy array of final weights of nodes
    """
    invD = np.array(invD)
    cA = np.array(A.toarray())


    W = np.matmul(cA, invD)


    # RWR
    alp0 = np.array(alpha * p0)
    p = alp0
    for i in range(iter):
        p = alp0 + (1 - alpha) * np.dot(W, p)


    return p

    # Attempt to do rwr using a present cupy library if available. If not, fall back to the numpy implementation
    # try:
    #     import cupy as cp
    #
    #     invD = cp.array(invD)
    #     cA = cp.array(A.toarray())
    #
    #     W = cp.matmul(cA, invD)
    #
    #     # RWR
    #     alp0 = cp.array(alpha * p0)
    #     p = alp0
    #     for i in range(iter):
    #         p = alp0 + (1 - alpha) * cp.dot(W, p)
    #
    #     return p
    #
    # except:


def graph_with_weights(wgraph, p, outcode='outputgraph', scale=True):
    """
    assigns its weight from the network propagation output (p) to each node.optionally will also scale the weights and
    skip labeling and adding weights to nodes that do not have a weight score higher than the cutoff.Writes the subgraph
    from the selected nodes into a gexf file and returns a string with the relative path of the written graph file

    :param wgraph: input graph
    :param p: weights to be added to the nodes
    :param outcode: name of the output gexf file
    :param scale: whether to scale the weights. Default = True
    :return: graph with weights added
    """

    # Rescale the weights between 0 and the scale_limit value
    if scale:   p *= scale_limit / p.max()

    sub_nodes = []
    i = 0
    for node in wgraph.nodes:  # for each node
        # add it's weight
        wgraph.nodes[node]['propagated_weight'] = p[i]

        i += 1


    return wgraph


def netprop(graph, seedlist, weight, alpha, iter, scale=True, regen=False):
    """
    Does network propagation in the given graph for the given seeds. The method can run for multiple times based on the
    number of values in the weight, alpha and iter lists (once for each combination). Optionally the weight values can
    be scaled and a cutoff value put for weights. The regen bool determines whether or not to regenerate the numpy
    arrays corresponding to a matrix used in the calculation. This can be set to false to allow matrices that do not
    change (ie degree matrix inverse, adjacency matrix etc) to be reused. This significantly reduces processing time
    when multiple values are used for weight/alpha/iter

    :param graph: input graph
    :param seedlist: path to seed file
    :param weight: list of default weight values to give when initializing network propagation
    :param alpha: list of alpha values (learning parameter) to give when initializing network propagation
    :param iter: list of iteration numbers to give when initializing network propagation
    :param scale: Whether to scale the value or not
    :param cutoff: defines the optional cutoff below which nodes will not be added to the output graph
    :param regen: regenerate the numpy arrays for some matrices on subsequent runs. Can usually be set to False
    :return: list of strings that are paths to created graphs
    """


    # Get adjacency matrix of graph
    A = nx.graphmatrix.adjacency_matrix(graph)


    # Create the inverse degree matrix and save as a temp file. If a file already exists and regen = False
    # then just use that instead
    try:
        if regen: raise FileNotFoundError  # if regen is set to True, force creation of file
        invD = np.load('tmp.npy')
    except FileNotFoundError:
        # Create degree matrix of graph
        D = np.eye(graph.number_of_nodes())

        i = 0
        for name, degree in graph.degree:
            D[i, i] = degree
            i += 1
        invD = np.linalg.inv(D)


    # Input weights for network
    p0 = weights_from_seeds(graph, seedlist, weight)
    # Do a random walk for some iterations and get the final weight vector for nodes
    p = _rwr(p0, alpha, A, invD, iter)

    outcode = f"{outfile}-w={weight}-a={alpha}-i={iter}"
    graphname = graph_with_weights(graph, p, outcode, scale)


    return graphname
