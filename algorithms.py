import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

def predict (ngraph,seeds,d,diff):

    def functional_flow(graph, seedlist, d):
        deg_matrix = csr_matrix((len(graph.nodes), 1)).toarray()  # array storing node degree
        rem_fluid = csr_matrix((len(graph.nodes), 1)).toarray()  # array storing remaining fluid in each reservoir
        en_fluid = csr_matrix((len(graph.nodes), 1)).toarray()  # array storing entered fluid in each node

        seeds = seedlist  # seeds in the network

        i = 0
        for name in graph.nodes:
            deg_matrix[i, 0] = graph.degree(name, "weight")  # preparing degree array
            if name in seeds:  # updating remaining fluid at t=0 time step
                rem_fluid[i, 0] = float("inf")  # if node is a seed infinite fluid exist
            i += 1

        adj_m = nx.adjacency_matrix(graph)  # getting adjacency matrix of the graph
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

            en_fluid = np.add(en_fluid, entered_fluid)  # update total enetered
            # print("rem", rem_fluid)

        n = 0

        for node in graph.nodes:  # assigning func score to each node
            graph.nodes[node]['functional_score'] = en_fluid[n, 0]
            n += 1

        return graph


    graphf = functional_flow(ngraph, seeds, d) # calculating functional score

    # to score each type of scores separately
    rwr={}
    fun = {}
    mv = {}
    hishigaki = {}

    scores = pd.DataFrame()  # to store all scores of predictions

    import network_prop  # calculating rwr


    graphr = network_prop.netprop(ngraph,seeds,100,0.1,5)

    for node in ngraph.nodes:  # for each node
        present = 0

        if node not in seeds:  # for unannotated nodes

            if node in diff:  # if validated by DEPs
                present += 1

            rwr[node] = graphr.nodes[node]['propagated_weight']
            fun[node] = graphf.nodes[node]['functional_score']

            li = [m for m in ngraph.neighbors(node)]  # getting list of neighbors
            known = [x for x in li if x in seeds]  # list of seeds as neighbors
            ef = (len(seeds) * len(li) / len(ngraph.nodes))
            nf = len(known)
            if nf != 0:
                score = (nf - ef) ** 2 / ef
                hishigaki[node] = score  # hishigaki score
            else:
                hishigaki[node] = 0
            mv[node] = len(known)  # mv score

            # tot = mv[node]+ hishigaki[node] + fun[node] + rwr[node] + present

            row = {"Node": node, "Majority voting score": mv[node], "Hishigaki score": hishigaki[node],
                   "Functional flow score": fun[node], "RWR": rwr[node], "Validated by DEPs": present}

            scores = scores.append(row, ignore_index=True)

        # for score normalization

        minrwr = min(scores["RWR"])
        maxrwr = max(scores["RWR"])
        minfun = min(scores["Functional flow score"])
        maxfun = max(scores["Functional flow score"])
        minmv = min(scores["Majority voting score"])
        maxmv = max(scores["Majority voting score"])
        minh = min(scores["Hishigaki score"])
        maxh = max(scores["Hishigaki score"])

        for index, row in scores.iterrows():
            rwr = (row["RWR"] - minrwr) / (maxrwr - minrwr)
            funsc = (row["Functional flow score"] - minfun) / (maxfun - minfun)
            mvsc = (row["Majority voting score"] - minmv) / (maxmv - minmv)
            hsc = (row["Hishigaki score"] - minh) / (maxh - minh)
            tot = rwr + funsc + mvsc + hsc + row["Validated by DEPs"]

            scores.at[index, "Normalized majority voting score"] = mvsc
            scores.at[index, "Normalized hishigaki score"] = hsc
            scores.at[index, "Normalized functional flow score"] = funsc
            scores.at[index, "Normalized rwr score"] = rwr
            scores.at[index, "Total Score"] = tot

        sorted_scores = scores.sort_values("Total Score", ascending=False)

        sorted_scores.to_excel("scores for predictions.xlsx")

















