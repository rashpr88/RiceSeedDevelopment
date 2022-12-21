import pandas as pd
import gzip

aliases = pd.read_table(gzip.open("4530.protein.aliases.v11.5.txt.gz", mode='rb'))

info = pd.read_table(gzip.open("4530.protein.info.v11.5.txt.gz"))

map_info = {}  # to store mapped preferred names
map_aliases = {}
data = pd.read_excel("seeds.xlsx", sheet_name="seeds")  # reading reference file as a dataframe
seeds = [m for m in data["preferredName"]]  # list of seeds

for i, row in info.iterrows():
	map_info[row['preferred_name']] = row['#string_protein_id']  # mapping b/w preferred name and string id
for i, row in aliases.iterrows():  # mapping b/w string id and ensemble gene id
	if aliases.at[i, 'source'] == 'Ensembl_gene':
		map_aliases[row['#string_protein_id']] = row['alias']

with open ("25/partitions.txt",'r') as file:
	all = ""
	for row in file:
		if row != "\n":
			string = ""
			index = row[:row.find("[") - 2]
			print(index)
			r = row[row.find("[") + 1:row.find("]")]
			for i in list(r.split(',')):
				if i != '':
					i = i[1:].replace("'", "")
					string = string + "\n" + map_aliases[map_info[i]]
					if i not in seeds:
						all = all + "\n" + map_aliases[map_info[i]]
			with open("25/" + str(index) + "-module.txt", "w") as wf:  # mapping ensemble gene ids in each module
				wf.write(string)
	with open("25/all.txt","w") as allm : # list of mapped candidates
		allm.write(all)
string = ""
with open ("25/hubs.txt",'r') as file:
	for row in file:
		row = row.strip('\n')
		if row != "" :
			string = string + "\n"+map_aliases[map_info[row]]
with open ("25/all_hubs.txt",'w') as file:  # mapping ensemble gene ids for hubs
	file.write(string)
