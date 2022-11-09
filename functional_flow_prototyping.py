import networkx as nwx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


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



#
# pos = {"s1":[4,6],"p1":[4,8],"p5":[2,6],"p4":[2,3],"p3":[6,3],"p2":[6,6],"e1":[8,6],"e2":[8,4],"s2":[10,6],"e3":[11,7],"e4":[12,8]}  # positions of nodes in the graph
# weight = nwx.get_edge_attributes(network,"weight")  # weights in graph
#
# col = ["yellow" if node in proteins else "blue" for node in network.nodes]  # seeds in yellow and rest in blue
#
# nwx.draw_networkx(network,pos=pos,with_labels=True,node_color = col)  # drawing the network
# nwx.draw_networkx_edge_labels(network,pos=pos,edge_labels=weight)  # adding edge labels
# plt.show()  # visualizing the graph drawn
# print(nwx.to_pandas_adjacency(network))
unknown = []  # unannotated proteins to the target function
for i in network.nodes: # filtering unannotated proteins to the target function
	if i not in proteins:
		unknown.append(i)

entered_fluid = {}  # to store fluid entered
remainder = {}  # to store remaining fluid
update_remainder = {}  # to update remainder at each node at the end of each time step

for p in unknown:
	entered_fluid[p] = 0
	remainder[p] = 0
	update_remainder[p] = 0

known_in = []
d = round(nwx.diameter(network)/2)
print(d,"d")

for i in proteins:
	if i in network.nodes:
		known_in.append(i)

n = [known_in]  # a nested list of nodes involving in flow
prev = []
a = nwx.to_numpy_array(network)
adde = False
addp = False

def update_fscore(e,p):
	print("score : ", e,"-", p)
	global adde
	global addp

	edge_w = network[e][p]["weight"]  # edge weight
	deg_p = network.degree(p, "weight")  # weighted degree for p
	deg_e = network.degree(e, "weight")  # weighted degree for e


	if p in known_in:  # for direct neighbors of seeds
		score = min(edge_w, (float("inf") * (edge_w / deg_p)))

		entered_fluid[e] = entered_fluid[e] + score
		update_remainder[e] = update_remainder[e] +score
		adde = True

	elif e in known_in:
		score = min(edge_w, (float("inf") * (edge_w / deg_e)))

		entered_fluid[p] = entered_fluid[p] + score
		update_remainder[p] = update_remainder[p] + score
		adde = True


	elif p not in known_in and e not in known_in:  # downhill flow
		if remainder[p] > remainder[e]:
			score = min(edge_w, (remainder[p] * (edge_w / deg_p)))
			update_remainder[p] = update_remainder[p] - score  # to update remainder in p
			entered_fluid[e] = entered_fluid[e] + score  # to update entered fluid in p
			update_remainder[e] = update_remainder[e] + score
			adde = True


		elif remainder[e] > remainder[p]:

			score = min(edge_w, (remainder[e] * (edge_w / deg_e)))
			update_remainder[e] = update_remainder[e] - score
			entered_fluid[p] = entered_fluid[p] + score
			update_remainder[p] = update_remainder[p] + score
			addp = True

	print("updating remaining volume in node p2",update_remainder["p2"])



for i in range(0, d):  # iterations based on graph radius
	interactions = []  # to track already checked interactions

	for s in range(0, len(n)):  # for each given set of nodes in each level

		new = []  # to store next level neighbors

		for p in n[s]:  # for proteins in each level involving the flow

			neigh = network.neighbors(p)  # getting list of neighbors of p protein

			for e in neigh:  # for each neighbor of p protein
				adde = False
				addp = False

				if (e in n[s-1] and s !=0) or ([p,e] in interactions or [e,p] in interactions) or (e in known_in and p in known_in):  # discard interactions already tracked in a previous level within same time step
					continue
				interactions.append([p, e])  # update checked interactions except for very 1st iteration

				update_fscore(e, p)  # to calculate and update functional score

				if s == len(n) - 1:  # updating n+1th level neighbors at nth level
					if adde is True:
						new.append(e)
					elif addp is True:
						new.append(p)


		if s == len(n) - 1:  # extending neighbor levels
			new = [*set(new)]  # eliminating duplicates
			n.append(new)
		print(n)

	for k in update_remainder.keys():  # updating remainder in neighbors when all interactions under iteration are covered
		remainder[k] = update_remainder[k]
	print("p2 entered fluid : ",entered_fluid["p2"])
	print("e12 entered fluid : ",entered_fluid["e1"])
	print("p2 remaining fluid : ",remainder["p2"])


print("Entered",entered_fluid)
print("Remainder",remainder)

sortd = sorted(entered_fluid.items(), key=lambda item: item[1], reverse=1)  # sorted dictionary
p = sortd[0]  # getting the protein record with highest functional flow score
result = "Protein\t\tScore\n"
for i in sortd:  # writing proteins and scores into a file
	result = result + str(i[0]) + "\t\t" + str(i[1]) + "\n"
f = open("Functional flow scorep.txt", "w")
f.write(result)
f.close()
print("Unknown protein with highest Functional flow score: ", str(p[0]) + " - " + str(p[1]))

