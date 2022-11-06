import pandas as pd
import networkx as nwx
import numpy as np


cand = set()

def predict(net):
	print("running algorithms")
	data = pd.read_excel("seeds.xlsx", sheet_name="seeds")  # reading reference file as a dataframe
	seeds = [m for m in data["preferredName"]]  # list of seeds

	nodes = [n for n in net.nodes]  # list of all nodes in the graph
	unknown = [n for n in nodes if n not in seeds]  # list of unknown proteins

	entered_fluid = {}  # to store fluid entered
	remainder = {}  # to store remaining fluid
	update_remainder = {}  # to update remainder in unannotated nodes

	for p in unknown:
		entered_fluid[p] = 0
		remainder[p] = 0
		update_remainder[p] = 0

	known_in = [s for s in seeds if s in nodes]  # known seeds in the network
	print(len(known_in))
	# d = nwx.radius(net) # half diameter of network
	# print("Number of iterations ",d)

	n = [known_in]  # a nested list of nodes involving in flow

	add = False  # decide whether to add neighbor for next level fluid flow

	def update_fscore(e, p):

		global add

		edge_w = net[e][p]["weight"]  # edge weight
		deg_p = net.degree(p, "weight")  # weighted degree for p
		deg_e = net.degree(e, "weight")  # weighted degree for e

		if p in known_in:  # for direct neighbors of seeds
			score = min(edge_w, (float("inf") * (edge_w / deg_p)))

			entered_fluid[e] = entered_fluid[e] + score
			update_remainder[e] = update_remainder[e] + score
			add = True

		elif e in known_in:
			score = min(edge_w, (float("inf") * (edge_w / deg_e)))

			entered_fluid[p] = entered_fluid[p] + score
			update_remainder[p] = update_remainder[p] + score


		elif p not in known_in and e not in known_in:  # downhill flow
			if remainder[p] > remainder[e]:
				score = min(edge_w, (remainder[p] * (edge_w / deg_p)))
				update_remainder[p] = update_remainder[p] - score  # to update remainder in p
				entered_fluid[e] = entered_fluid[e] + score  # to update entered fluid in p
				update_remainder[e] = update_remainder[e] + score
				add = True


			elif remainder[e] > remainder[p]:

				score = min(edge_w, (remainder[e] * (edge_w / deg_e)))
				update_remainder[e] = update_remainder[e] - score
				entered_fluid[p] = entered_fluid[p] + score
				update_remainder[p] = update_remainder[p] + score


	for i in range(0, 5):  # iterations based on graph radius
		interactions = []  # to track already checked interactions

		for s in range(0, len(n)):  # for each given set of nodes in each level

			new = []  # to store next level neighbors

			for p in n[s]:  # for proteins in each level involving the flow

				neigh = net.neighbors(p)  # getting list of neighbors of p protein

				for e in neigh:  # for each neighbor of p protein
					add = False

					if (e in n[s - 1] and s != 0) or ([p, e] in interactions or [e, p] in interactions) or (
							e in known_in and p in known_in):  # discard interactions already tracked in a previous level within same time step
						continue
					interactions.append([p, e])  # update checked interactions except for very 1st iteration

					update_fscore(e, p)  # to calculate and update functional score

					if s == len(n) - 1 and add is True:  # updating n+1th level neighbors at nth level
						new.append(e)

			if s == len(n) - 1:  # extending neighbor levels
				new = [*set(new)]  # eliminating duplicates
				n.append(new)

		for k in update_remainder.keys():  # updating remainder in neighbors when all interactions under iteration are covered
			remainder[k] = update_remainder[k]

	print("Entered", entered_fluid)
	print("Remainder", remainder)

	votes = {}  # to store proteins and corresponding majority scores
	hishigaki = {}

	for i in unknown:  # calculating hishigaki and majority voting scores
		li = [m for m in net.neighbors(i)]  # getting list of neighbors
		known = [x for x in li if x.upper() in known_in]  # list of seeds as neighbors
		ef = (len(known_in) * len(li) / len(net.nodes))  # expected frequency for the function
		nf = len(known)
		if nf != 0:
			score = (nf - ef) ** 2 / ef
			hishigaki[i] = score
		votes[i] = len(known) # mv score



	def splitter(n):  # to get common candidates above 80th percentile
		global cand
		candidates = {}
		for i in range(0, len(n)):
			candidates[n[i][0]] = n[i][1]
		array = list(candidates.values())
		boundary = np.percentile(array, 80)

		s = set()

		for i in candidates.keys():  # filtering candidates in 80th percentile
			if float(candidates[i]) > boundary:
				s.add(i)
		if len(cand) == 0:
			cand = s

		else:
			cand = cand.intersection(s)

	def sorting(s1, n1, s2, n2, s3, n3):

		for i, name in zip((s1, s2, s3), (n1, n2, n3)):
			sortd = sorted(i.items(), key=lambda item: item[1], reverse=1)  # sorted dictionary
			p = sortd[0]  # getting the protein record with highest functional flow score


			result = "Protein\t\tScore\n"

			for i in sortd:  # writing proteins and scores into a file
				result = result + str(i[0]) + "\t\t" + str(i[1]) + "\n"
			f = open(name + ".txt", "w")
			f.write(result)
			f.close()
			print("Unknown protein with highest " + name + " : ", str(p[0]) + " - " + str(p[1]))

			splitter(sortd)

	sorting(entered_fluid, "Functional flow score", hishigaki, "Hishigaki score", votes, "Majority voting score")

	print("candidates supported by 3 algorithms : ", cand)
	cand_res = "Candidates supported by 3 algorithms\n\n"
	for i in cand:
		cand_res = cand_res + "\n" + str(i)

	f2 = open("candidates.txt", "w")
	f2.write(cand_res)
	f2.close()

	nodes_for_subnetwork = cand.union((set(known_in)))  # extracting rice seed development sub-network

	for edge in net.edges:
		if not nodes_for_subnetwork.intersection(set(edge)) == set(edge):
			net.remove_edge(edge[0], edge[1])

	nwx.write_gml(net, "rice_seed_development_sub_network.gml")  # sub-network extracted




