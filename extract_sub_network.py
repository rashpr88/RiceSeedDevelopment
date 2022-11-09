import os
import networkx as nwx
import pandas as pd


if os.path.exists("./rice_network.gml"):  # if the created network already exists
	print(True)
	networkp = nwx.read_gml("rice_network.gml")

	if os.path.exists("scores for predictions.xlsx"):
		data = pd.read_excel("scores for predictions.xlsx",sheet_name="Sheet1")
		candidates = [c for c in data["Node"]]

		candidates = candidates[0:100]  # top 100 candidates
		for edge in networkp.edges:
			if not set(candidates).intersection(set(edge)) == set(edge):
				networkp.remove_edge(edge[0], edge[1])

		nwx.write_gml(networkp, "rice_seed_development_sub_network.gml")  # sub-network extracted








