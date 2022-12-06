import networkx as nwx
import pandas as pd
import gzip
import os

networkp = nwx.Graph()  # initiating the network graph
data = pd.read_excel("seeds.xlsx", sheet_name="seeds")  # reading reference file as a dataframe
seeds = [m for m in data["preferredName"]]  # list of seeds

if os.path.exists("./rice_network.gml"):  # if the created network already exists
	rice_net = nwx.read_gml("rice_network.gml")
else:
	map = {}  # to store mapped preferred names

	info = pd.read_table(gzip.open("4530.protein.info.v11.5.txt.gz", mode='rb'))

	for i, row in info.iterrows():
		map[row['#string_protein_id']] = row['preferred_name']

	data = pd.read_table(gzip.open("4530.protein.links.v11.5.txt.gz", mode='rb'), sep=" ")
	data = data.drop_duplicates()  # reading interaction file as a dataframe

	for index, row in data.iterrows():  # replacing string ids with preferred names
		p1 = row["protein1"]
		p2 = row["protein2"]
		pr1 = map[p1]
		pr2 = map[p2]
		sc = row["combined_score"] / 1000

		if sc >= 0.4:  # to filter interactions based on combined score
			networkp.add_edge(pr1, pr2, weight=sc)

	rice_net = nwx.Graph()

	if not nwx.is_connected(networkp):
		sub_graphs = max(nwx.connected_components(networkp), key=len)
		rice_net = networkp.subgraph(sub_graphs)
	else:
		rice_net = networkp

	for node in rice_net.nodes:  # for seed visualization
		if node in seeds:
			rice_net.nodes[node]["seeds"] = 1
		else:
			rice_net.nodes[node]["seeds"] = 0


	nwx.write_gml(rice_net, "rice_network.gml")  # writing prepared network to a gml file

known_in = [s for s in seeds if s in rice_net.nodes]  # known seeds in the network
print("No of seeds in the network : ",len(known_in))
# d = nwx.radius(rice_net) # half diameter of network
d = 5
# print("Number of iterations : ",d)

deps = pd.read_excel("seeds.xlsx", sheet_name="DEPs")
diff = [m for m in deps["preferredName"]]  # list of DEPs

import algorithms

algorithms.predict(rice_net,known_in,d,diff)  # to predict candidates

