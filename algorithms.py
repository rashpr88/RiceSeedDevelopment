import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse

def predict (ngraph,seeds,d,diff):

    def functional_flow(graph, seedlist, d):
        deg_matrix = np.zeros((len(graph.nodes), 1))  # array storing node degree
        rem_fluid = np.zeros((len(graph.nodes), 1))  # array storing remaining fluid in each reservoir
        en_fluid = np.zeros((len(graph.nodes), 1))  # array storing entered fluid in each node

        i = 0
        for name in graph.nodes:
            deg_matrix[i, 0] = 1 / (graph.degree(name, "weight"))  # preparing degree array
            if name in seedlist:  # updating remaining fluid at t=0 time step
                rem_fluid[i, 0] = float("inf")  # if node is a seed infinite fluid exist
                en_fluid[i, 0] = float("inf")  # if node is a seed infinite fluid pumped in
            i += 1

        adj_matrix = nx.adjacency_matrix(graph)  # getting adjacency matrix of the graph
        adj_matrix = sparse.csr_matrix(adj_matrix.toarray())  # converting to a sparse matrix

        for i in range(0, d):

            t_rem_fluid = rem_fluid.transpose()  # transpose of remaining fluid

            mesh_rem_fluid = np.meshgrid(rem_fluid, t_rem_fluid)  # comparing remainder in each reservoir

            max_rem_fluid = np.maximum(*mesh_rem_fluid)  # picking reservoir with max remainder

            greater = np.greater(t_rem_fluid, rem_fluid)  # getting a boolean array for remainder comparison
            greater = sparse.csr_matrix(greater * 1)  # to track whether a replacement happens

            less = np.less(t_rem_fluid, rem_fluid)  # to distinguish between equal and less
            less = sparse.csr_matrix(less * -1)

            positives = sparse.csr_matrix.multiply(greater,adj_matrix)

            negatives = sparse.csr_matrix.multiply(less, adj_matrix)

            n = negatives.copy()

            p = positives.copy()

            n[n < 0] = 1  # to update deg matrix for fluid entering interactions

            p[p > 0] = 1  # to update deg matrix for fluid leaving interactions

            res = []

            def score_calculation(sign1, sign_list1, deg1, clue1, sign2, sign_list2, deg2, clue2):

                for s, l, d, t in zip((sign1, sign2), (sign_list1, sign_list2), (deg1, deg2), (clue1, clue2)):

                    updated_deg_matrix = sparse.csr_matrix.multiply(s,d)  # degree matrix for the interactions under consideration

                    weight_proportion = sparse.csr_matrix.multiply(l,updated_deg_matrix)  # expressing flow capacity as a proportion

                    reservoir_cap = sparse.csr_matrix.multiply(s,max_rem_fluid)  # updating current fluid volumes at reservoirs

                    using_remainder_and_weights = sparse.csr_matrix.multiply(weight_proportion,reservoir_cap)  # calculating possible flow volume

                    if t == "p":
                        fluid_v = sparse.csr_matrix.minimum(l,using_remainder_and_weights)  # picking the minimum out of the edge degree and calculated flow volume

                    elif t == "n":

                        fluid_v= sparse.csr_matrix.maximum(l,using_remainder_and_weights)  # picking the minimum out of the edge degree and calculated flow volume

                    res.append(fluid_v)  # adding to result

            score_calculation(p, positives, deg_matrix.transpose(), "p", n, negatives, deg_matrix,
                              "n")  # calculating functional score

            entered_fluid = sparse.csr_matrix.sum(res[0], axis=1)  # entered fluid in ith iteration

            fluid_exit = sparse.csr_matrix.sum(res[1], axis=1)  # fluid exit in ith iteration

            update_rem = np.add(entered_fluid, fluid_exit)  # remainder due to ith iteration

            rem_fluid = np.add(rem_fluid, update_rem)  # update total remainder

            en_fluid = np.add(en_fluid, entered_fluid)  # update total entered

        x = 0

        for node in graph.nodes:  # assigning func score to each node
            graph.nodes[node]['functional_score'] = en_fluid[n, 0]
            x += 1

        return graph


    functional_flow(ngraph, seeds, d) # calculating functional score

    scores = pd.DataFrame()  # to store all scores of predictions

    import network_prop  # calculating rwr

    network_prop.netprop(ngraph,seeds,100,0.1,5)

    for node in ngraph.nodes:  # for each node
        present = 0

        if node not in seeds:  # for unannotated nodes

            if node in diff:  # if validated by DEPs
                present += 1

            rwr = ngraph.nodes[node]['propagated_weight']
            fun = ngraph.nodes[node]['functional_score']

            li = [m for m in ngraph.neighbors(node)]  # getting list of neighbors
            known = [x for x in li if x in seeds]  # list of seeds as neighbors
            ef = (len(seeds) * len(li) / len(ngraph.nodes))
            nf = len(known)
            if nf != 0:
                score = (nf - ef) ** 2 / ef
                hishigaki = score  # hishigaki score
            else:
                hishigaki = 0
            mv = len(known)  # mv score

            row = {"Node": node, "Majority voting score": mv, "Hishigaki score": hishigaki,
                   "Functional flow score": fun, "RWR": rwr, "Validated by DEPs": present}

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
        tot = rwr + funsc + mvsc + hsc

        scores.at[index, "Normalized majority voting score"] = mvsc
        scores.at[index, "Normalized hishigaki score"] = hsc
        scores.at[index, "Normalized functional flow score"] = funsc
        scores.at[index, "Normalized rwr score"] = rwr
        scores.at[index, "Total Score"] = tot

    sorted_scores = scores.sort_values("Total Score", ascending=False)

    sorted_scores.to_excel("scores for predictions.xlsx")


















