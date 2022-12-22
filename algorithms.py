import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse


def _array_initiation_(graph, seeds):  # preparing required arrays
	degree_matrix = np.zeros((len(graph.nodes), 1))  # array storing node degree
	remaining_fluid = np.zeros((len(graph.nodes), 1))  # array storing remaining fluid in each reservoir
	entered_fluid = np.zeros((len(graph.nodes), 1))  # array storing entered fluid in each node

	i = 0
	for name in graph.nodes:
		degree_matrix[i, 0] = 1 / (graph.degree(name, "weight"))  # preparing degree array
		if name in seeds:  # updating remaining fluid at t=0 time step
			remaining_fluid[i, 0] = float("inf")  # if node is a seed infinite fluid exist
			entered_fluid[i, 0] = float("inf")  # if node is a seed infinite fluid pumped in
		i += 1

	adj_matrix = nx.adjacency_matrix(graph)  # getting adjacency matrix of the graph
	adj_matrix = sparse.csr_matrix(adj_matrix.toarray())  # converting to a sparse matrix

	return [adj_matrix, degree_matrix, entered_fluid, remaining_fluid]


def score_calculation(sign1, sign_list1, deg1, clue1, sign2, sign_list2, deg2, clue2, max_rem_fluid):  # calculating functional flow scores
	res = []  # to store updated scores as two separate arrays for fluid entering and leaving interactions
	for s, l, d, t in zip((sign1, sign2), (sign_list1, sign_list2), (deg1, deg2), (clue1, clue2)):

		updated_deg_matrix = sparse.csr_matrix.multiply(s, d)  # degree matrix for the interactions under consideration

		weight_proportion = sparse.csr_matrix.multiply(l, updated_deg_matrix)  # flow capacity as a proportion of weights

		reservoir_cap = sparse.csr_matrix.multiply(s, max_rem_fluid)  # updating current fluid volumes at reservoirs

		using_remainder_and_weights = sparse.csr_matrix.multiply(weight_proportion,reservoir_cap)  # calculating the volume allowed to flow

		if t == "p":
			fluid_volume = sparse.csr_matrix.minimum(l,using_remainder_and_weights)  # picking the minimum out of the edge degree and allowed volume to flow

		elif t == "n":

			fluid_volume = sparse.csr_matrix.maximum(l,using_remainder_and_weights)  # picking the minimum out of the edge degree and allowed volume to flow

		res.append(fluid_volume)  # adding to result
	return res


def functional_flow(graph, seeds, d):
	arrays = _array_initiation_(graph, seeds)

	for i in range(0, d):
		t_rem_fluid = arrays[3].transpose()  # transpose of remaining fluid

		mesh_rem_fluid = np.meshgrid(arrays[3], t_rem_fluid)  # comparing remainder in each reservoir

		max_rem_fluid = np.maximum(*mesh_rem_fluid)  # picking reservoir with maximum fluid remaining

		greater = np.greater(t_rem_fluid, arrays[3])  # getting a boolean array for remainder comparison
		greater = sparse.csr_matrix(greater * 1)  # to track maximum values

		less = np.less(t_rem_fluid, arrays[3])  # to distinguish between equal and less
		less = sparse.csr_matrix(less * -1)

		positives = sparse.csr_matrix.multiply(greater, arrays[0])  # adjacency matrix for fluid entering reactions

		negatives = sparse.csr_matrix.multiply(less, arrays[0])  # adjacency matrix for fluid leaving interactions

		n = negatives.copy()

		p = positives.copy()

		n[n < 0] = 1  # to update deg matrix for fluid entering interactions

		p[p > 0] = 1  # to update deg matrix for fluid leaving interactions

		result = score_calculation(p, positives, arrays[1].transpose(), "p", n, negatives, arrays[1],"n", max_rem_fluid)

		entered_fluid = sparse.csr_matrix.sum(result[0], axis=1)  # entered fluid in i(th) iteration

		fluid_exit = sparse.csr_matrix.sum(result[1], axis=1)  # fluid exit during i(th) iteration

		update_rem = np.add(entered_fluid, fluid_exit)  # remaining fluid due to ith iteration

		arrays[3] = np.add(arrays[3], update_rem)  # update total remaining fluid

		arrays[2] = np.add(arrays[2], entered_fluid)  # update total entered fluid

	x = 0

	for node in graph.nodes:  # assigning calculated functional flow score to each node
		graph.nodes[node]['functional_score'] = arrays[2][x, 0]
		x += 1

	return graph


def predict(graph, seeds, d, diff):  # getting the ensemble prediction
	functional_flow(graph, seeds, d)  # calculating functional flow score

	import network_prop

	network_prop.netprop(graph, seeds, 100, 0.1, 5)  # calculating rwr

	scores = pd.DataFrame()  # to store all scores of predictions

	for node in graph.nodes:  # for each node
		present = 0

		if node not in seeds:  # for unannotated nodes

			if node in diff:  # if validated by DEPs
				present += 1

			rwr = graph.nodes[node]['propagated_weight']
			fun = graph.nodes[node]['functional_score']

			li = [m for m in graph.neighbors(node)]  # getting list of neighbors
			known = [x for x in li if x in seeds]  # list of seeds as neighbors
			ef = (len(seeds) * len(li) / len(graph.nodes))
			nf = len(known)
			if nf != 0:
				score = (nf - ef) ** 2 / ef
				hishigaki = score  # hishigaki score
			else:
				hishigaki = 0
			mv = len(known)  # mv score

			row = {"Node": node, "Majority voting score": mv, "Hishigaki score": hishigaki,"Functional flow score": fun, "RWR": rwr, "Validated by DEPs": present}

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
		tot = rwr + funsc + mvsc + hsc + scores.at[index, "Validated by DEPs"]

		scores.at[index, "Normalized majority voting score"] = mvsc
		scores.at[index, "Normalized hishigaki score"] = hsc
		scores.at[index, "Normalized functional flow score"] = funsc
		scores.at[index, "Normalized rwr score"] = rwr
		scores.at[index, "Total Score"] = tot

	sorted_scores = scores.sort_values("Total Score", ascending=False)

	sorted_scores.to_excel("scores for predictions_0.7.xlsx")


















