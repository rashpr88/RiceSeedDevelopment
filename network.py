import networkx as nwx
import pandas as pd
import gzip
import os

networkp = nwx.Graph()  # initiating the network graph

if os.path.exists("./rice_network.gml"):  # if the created network already exists
	print(True)
	networkp = nwx.read_gml("rice_network.gml")
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


print("no of edges network before pruning",networkp.number_of_edges())
prune_p = nwx.Graph()

if not nwx.is_connected(networkp) :
	print("not connected")
	sub_graphs = max(nwx.connected_components(networkp),key=len)
	prune_p = networkp.subgraph(sub_graphs)
else:
	prune_p = networkp

nwx.write_gml(prune_p, "rice_network.gml")  # writing prepared network to a gml file
print("no of edges network after pruning", prune_p.number_of_edges())

import algorithms

algorithms.predict(prune_p)  # to predict candidates



# proteins = data.protein1.tolist() + data.protein2.tolist()  # list to search for preferred names
# proteins = list(set(proteins))  # without duplicates


# for i,row in info.iterrows():
# 	map[row['#string_protein_id']] = row['preferred_name']
# num = 4000 * (math.ceil(len(proteins) / 4000)) # to replace preferred name directly from db
#
# i = 4000
# n = 0
# x = 0
# while (i < num+1) :  # mapping preferred names with string ids
# 	if len(proteins) - (i-4000) > 4000:
# 		x = i
# 	else:
# 		x = (i-4000) + (len(proteins) - (i-4000))
# 	string_ids = stringdb.get_string_ids(proteins[n:x], species=4530)
# 	# string_ids = string_ids.drop_duplicates(subset="stringId", keep="first")  # eliminating duplicates
# 	for l, r in string_ids.iterrows():  # evidence and process from uniprot
# 		id = string_ids.at[l, 'stringId']
# 		preferred_name = string_ids.at[l, "preferredName"]
# 		map[id] = preferred_name
#
# 	n = x
#
# 	i = i+ 4000



# f = open("interactions_updated.txt","w")
# result = data.to_string(header=True)
# f.write(result)

