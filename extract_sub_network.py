import os
import networkx as nwx
import pandas as pd
from community import community_louvain
import numpy as np

if os.path.exists("./rice_network.gml"):  # if the created network already exists
	print(True)
	networkp = nwx.read_gml("rice_network.gml")
	subnet = nwx.Graph()

	if os.path.exists("scores for predictions.xlsx"):
		data = pd.read_excel("scores for predictions.xlsx",sheet_name="Sheet1")
		candidates = [c for c in data["Node"]]
		data_s = pd.read_excel("scores for predictions.xlsx",sheet_name="seeds")
		seeds = [s for s in data_s["preferredName"]]


		candidates = candidates[0:100]  # top 100 candidates
		nodes = seeds + candidates
		for edge in networkp.edges:
			if not set(nodes).intersection(set(edge)) == set(edge):
				subnet.add_edge(edge[0],edge[1],weight=networkp[edge[0]][edge[1]]["weight"])

		nwx.write_gml(subnet, "rice_seed_development_sub_network.gml")  # sub-network extracted


partitions = community_louvain.best_partition(subnet)
p = set(partitions.values())

for node in partitions:
	subnet.nodes[node]["cluster"] = partitions[node]

l = subnet.degree(subnet.nodes, "weight")
deg = [ i[1] for i in l]
boundary = np.percentile(deg, 90)

for node in subnet.nodes:
	degree = subnet.degree(node, "weight")
	inter= set()

	for edge in subnet.edges:
		if node in edge :
			if partitions[edge[0]] != partitions[edge[1]]:
				inter.update([partitions[edge[0]],partitions[edge[0]]])
	inter.discard(partitions[node])
	print(inter)
	print(partitions[node])
	if len(inter) >= p/2:
		if degree > boundary:
			subnet.nodes[node]["hub"] = 3
		else:
			subnet.nodes[node]["hub"] = 1
	else:
		if degree > boundary:
			subnet.nodes[node]["hub"] = 2
		else:
			subnet.nodes[node]["hub"] = 0

nwx.write_gml(subnet,"sub_modules.gml")






