import os
import networkx as nwx
import pandas as pd
import numpy as np
from community import community_louvain
from collections import defaultdict

if os.path.exists("./rice_network.gml"):  # if the created network already exists
	print(True)
	networkp = nwx.read_gml("rice_network.gml")
	subnet = nwx.read_gml("rice_network.gml")

	if os.path.exists("scores for predictions.xlsx"):
		data = pd.read_excel("scores for predictions.xlsx",sheet_name="Sheet1")
		candidates = [c for c in data["Node"]]
		data_s = pd.read_excel("seeds.xlsx", sheet_name="seeds")
		seeds = [s for s in data_s["preferredName"] if s in networkp.nodes]


		candidates = candidates[0:100]  # top 100 candidates
		nodes = seeds + candidates
		subnet = networkp.subgraph(nodes)

		partitions = community_louvain.best_partition(subnet)
		result = defaultdict(list)
		avg ={}
		std = {}
		intra_deg = {}
		string = ""

		for key, val in sorted(partitions.items()):  # grouping neighbors based on partition
			result[val].append(key)

		for key, val in sorted(result.items()):
			string = string + "\n " + str(key) + ' : ' + str(val)
			interactions = []
			for i in val :
				module_neighbors = [m for m in subnet.neighbors(i) if m in val]
				intra_deg[i] = len(module_neighbors)
				interactions.append(len(module_neighbors))
			print(interactions)


			ave = np.mean(interactions)
			print(ave)
			avg[key] = ave
			s_deviation = np.std(interactions)
			print(s_deviation)
			std [key] = s_deviation

		with open("partitions.txt",'w') as f:
			f.write(string)
		z_tot = []
		hubs = ""

		for node in partitions:  # getting partition coefficient for each node
			module = partitions[node]
			subnet.nodes[node]["cluster"] = module


			pc = 0
			if std[module] != 0:
				z = (intra_deg[node] - avg[module]) / std[module]
				z_tot.append(z)

			else :
				z = 0
				z_tot.append(z)

			if z >= 0.9 :
				degree = subnet.degree(node)  # degree of the node
				neighbors = [m for m in subnet.neighbors(node)]
				neighbor_partitions = {key: partitions[key] for key in neighbors}

				result = defaultdict(list)

				for key, val in sorted(neighbor_partitions.items()):  # grouping neighbors based on partition
					result[val].append(key)

				for key, val in result.items():  # calculating PC
					pc += (len(val) / degree) ** 2
				pc = 1 - pc

				if pc > 0.5 :
					pc = 2
				else :
					pc = 1
				hubs = hubs + "\n"+node+ "\n"+ str(result) +"\n module: " +str(module)+"\n modules connecting : "+str(result.keys())+ "\n" + str(pc)+ "\n\n"

			else:
				pc = 0
			subnet.nodes[node]["pc"] = pc  # adding the node PC value

		with open("hubs.txt", 'w') as f:
			f.write(hubs)


		boundary = np.percentile(z_tot,90)
		print(boundary)

		nwx.write_gml(subnet, "rice_seed_development_sub_network.gml")  # partitioned sub-network extracted
