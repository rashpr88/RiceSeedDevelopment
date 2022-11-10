import networkx as nwx
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
from scipy import sparse

network = nwx.Graph()

proteins = ["s1",'s2']
network.add_edge("s1","p1",weight=0.3)
network.add_edge("s1","p2",weight=0.2)
network.add_edge("s1","p3",weight=0.5)
network.add_edge("s1","p4",weight=0.7)
network.add_edge("s1","p5",weight=0.6)
network.add_edge("p2","e1",weight=0.1)
network.add_edge("p2","e2",weight=0.4)
network.add_edge("e1","s2",weight=0.65)
network.add_edge("s2","e3",weight=0.7)
network.add_edge("e3","e4",weight=0.8)


unknown = []  # unannotated proteins to the target function
for i in network.nodes: # filtering unannotated proteins to the target function
	if i not in proteins:
		unknown.append(i)

d = nwx.radius(network)
print(d)
"""
Attempt to code functional flow
"""


def functional_flow(graph, seedlist, d):
    deg_matrix = np.zeros((len(graph.nodes), 1))  # array storing node degree
    rem_fluid = np.zeros((len(graph.nodes), 1))  # array storing remaining fluid in each reservoir
    en_fluid = np.zeros((len(graph.nodes), 1))  # array storing entered fluid in each node

    seeds = seedlist  # seeds in the network

    i = 0
    for name in graph.nodes:
        deg_matrix[i, 0] = graph.degree(name, "weight")  # preparing degree array
        if name in seeds:  # updating remaining fluid at t=0 time step
            rem_fluid[i, 0] = float("inf")  # if node is a seed infinite fluid exist
        i += 1

    adj_m = nwx.adjacency_matrix(graph)  # getting adjacency matrix of the graph
    adj_matrix = np.array(adj_m.toarray())  # converting to an array

    for i in range(0, d):

        t_rem_fluid = rem_fluid.transpose()  # transpose of remaining fluid

        mesh_rem_fluid = np.meshgrid(rem_fluid, t_rem_fluid)  # comparing remainder in each reservoir

        max_rem_fluid = np.maximum(*mesh_rem_fluid)  # picking reservoir with max remainder

        if_replaced = np.greater(t_rem_fluid, rem_fluid)  # getting a boolean array for remainder comparison
        onezero = np.multiply(if_replaced, 1)  # to track whether a replacement happens
        less = np.less(t_rem_fluid, rem_fluid)  # to distinguish between equal and less
        less = np.multiply(less, -1)

        up = np.add(less, onezero)  # array indicating flow direction

        current_interactions = np.multiply(up, adj_matrix)  # interactions involving the flow in ith iteration
        current_interactions[current_interactions == -0] = 0  # omitting zeros with sign
        # current_interactions[current_interactions == 0] = float("nan")

        positives = current_interactions.copy()
        positives[positives < 0] = 0  # fluid entering interactions

        negatives = current_interactions.copy()
        negatives[negatives > 0] = 0  # fluid leaving interactions

        n = negatives.copy()
        p = positives.copy()

        n[n < 0] = 1  # to update deg matrix for fluid entering interactions
        p[p > 0] = 1  # to update deg matrix for fluid leaving interactions

        res = []

        def score_calculation(sign1, sign_list1, deg1, clue1, sign2, sign_list2, deg2, clue2):

            for s, l, d, t in zip((sign1, sign2), (sign_list1, sign_list2), (deg1, deg2), (clue1, clue2)):
                s[s == 0] = float("nan")
                l[l == 0] = float("nan")

                updated_deg_matrix = np.multiply(s, d)  # degree matrix for the interactions under consideration

                weight_proportion = np.divide(l, updated_deg_matrix)  # expressing flow capacity as a proportion

                reservoir_cap = np.multiply(max_rem_fluid, s)  # updating current fluid volumes at reservoirs

                using_remainder_and_weights = np.multiply(reservoir_cap,
                                                          weight_proportion)  # calculating possible flow volume
                # print(using_remainder_and_weights)

                if t == "p":
                    mat = np.minimum(l,
                                     using_remainder_and_weights)  # picking the minimum out of the edge degree and calculated flow volume
                elif t == "n":

                    mat = np.maximum(l,
                                     using_remainder_and_weights)  # picking the minimum out of the edge degree and calculated flow volume

                # print("m", mat)
                mat = np.nan_to_num(mat)  # replacing nan with 0

                res.append(mat)  # adding to result

        score_calculation(p, positives, deg_matrix.transpose(), "p", n, negatives, deg_matrix,
                          "n")  # calculating functional score

        # print("res",res)

        entered_fluid = np.matrix(res[0]).sum(axis=1)  # entered fluid in ith iteration
        # print("enp", entered_fluid)

        fluid_exit = np.matrix(res[1]).sum(axis=1)  # fluid exit in ith iteration
        # print("ex", fluid_exit)

        update_rem = np.add(entered_fluid, fluid_exit)  # remainder due to ith iteration

        rem_fluid = np.add(rem_fluid, update_rem)  # update total remainder

        en_fluid = np.add(en_fluid, entered_fluid)  # update total entered
        # print("rem", rem_fluid)

    n = 0

    for node in graph.nodes:  # assigning func score to each node
        graph.nodes[node]['functional_score'] = en_fluid[n, 0]
        n += 1

    return graph


fun ={}
g = functional_flow(network,proteins,d)


for node in g.nodes:  # for each node
	# add it's weight
	if node not in proteins:
		fun[node] = g.nodes[node]['functional_score']
print(fun)
sortd = sorted(fun.items(), key=lambda item: item[1], reverse=1)  # sorted dictionary
print(sortd)