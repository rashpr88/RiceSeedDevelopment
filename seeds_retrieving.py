import pandas as pd
import stringdb
from Bio import Entrez
import requests, json
import numpy as np


data = pd.read_excel("ref.xlsx","seed_collection") # reading seeds as a dataframe
data = data.drop_duplicates() # dropping duplicates
names = ("hypothetical protein","unnamed protein product","hypothetical","protein of unknown")  # unambiguous
data = data[~data["Protein Name"].isin(names)] # to remove unambiguous proteins in the list

Entrez.email = "Your.Name.Here@example.org"

for i in set(data["Protein Name"]):  # to get protein name using ncbi id provided
	id = data.loc[data["Protein Name"]==i]["ID"] # to extract id
	index = id.to_string().split(" ")[0]
	id = id.to_string().split(" ")[4]

	if id != "NaN"  : # to get exact protein name using given ncbi id
		try:
			handle = Entrez.esummary(db="protein", id=id)
			record = Entrez.read(handle)
			title = record[0]["Title"]
			name = title[:title.index("[")]
			data.at[index,"Protein Name"] = name  # replacing the name with name in ncbi
		except:
			continue
string_ids = stringdb.get_string_ids(data["Protein Name"],species=4530)  # getting stringdb info

string_ids = string_ids.drop_duplicates(subset="stringId", keep="first")  # eliminating duplicates

api = "https://rest.uniprot.org/"


for i,r in string_ids.iterrows():  # evidences and biological processes from uniprot
	entry = string_ids.at[i,'stringId']
	r = requests.get(f"{api}/uniprotkb/search?query="+entry+" AND (taxonomy_id:4530)&size=1")
	record = json.loads(r.text)  # obtaining as a dictionary
	x = record["results"]


	if x != [] :
		acc = x[0]["primaryAccession"]
		y = x[0]["uniProtKBCrossReferences"]
		evi = ""  # to store evidences other than IEA
		process = ""  # to store process
		evic =""  # store IEA
		processc = ""  # computationally annotated process

		for k in y:
			if k["database"] == "GO" and "P:" in k['properties'][0]["value"]:  # extracting processes and evidences
				if "IEA" not in k['properties'][1]["value"]:  # not IEA
					if evi != "":
						evi = evi + "," + str(k['properties'][1]["value"])
						process = process + "," + str(k['properties'][0]["value"][2:])
					else:
						evi = evi + str(k['properties'][1]["value"])
						process = process + str(k['properties'][0]["value"][2:])
				else:
					if evic != "": # computational annotations
						evic = evic + "," + str(k['properties'][1]["value"])
						processc = processc + "," + str(k['properties'][0]["value"][2:])
					else:
						evic = evic + str(k['properties'][1]["value"])
						processc = processc + str(k['properties'][0]["value"][2:])

		string_ids.at[i, 'Uniprot accession'] = acc  # adding uniprot accession
		if evi != "":  # when manual evidences exist
			string_ids.at[i, 'Evidences'] = evi
			string_ids.at[i, "Biological Processes"] = process
		else:
			if evic != "":  # when computational evidences exist
				string_ids.at[i, 'Evidences'] = evic
				string_ids.at[i, "Biological Processes"] = processc
			else:
				string_ids.at[i, 'Evidences'] = np.nan
				string_ids.at[i, "Biological Processes"] = np.nan

no_annotations = string_ids[string_ids["Evidences"].isna()]  # no uniprot evidences
computational = string_ids.dropna()[string_ids.dropna()["Evidences"].str.contains("IEA")]  # computational annotations
seeds = string_ids.dropna()[~string_ids.dropna()["Evidences"].str.contains("IEA")]  # no computational annotations

query2 = set(set(data["Protein Name"]) - set(string_ids["queryItem"]))  # proteins not in string
not_in_string = pd.DataFrame(data=query2, columns=["Protein Name"])

deps = pd.read_excel("ref.xlsx","DEPs") # reading DEPs as a dataframe
deps_string_ids = stringdb.get_string_ids(deps['Protein Name'],species=4530)  # getting stringdb info

deps_string_ids = deps_string_ids.drop_duplicates(subset="stringId", keep="first")  # eliminating duplicates

seeds_in_deps = set(deps_string_ids["stringId"]).intersection(set(seeds["stringId"]))
print(len(seeds_in_deps))

for i,row in deps_string_ids.iterrows():  # removing seeds which are DEPs
	if row["stringId"] in seeds_in_deps:
		deps_string_ids.drop(i,axis=0,inplace=True)



deps_string_ids = string_ids.drop_duplicates(subset="stringId", keep="first")  # eliminating duplicates
with pd.ExcelWriter ("seeds.xlsx") as w:  # writing into separate excel sheets
	seeds.to_excel(w, sheet_name="seeds")  # seeds with manual annotations
	computational.to_excel(w, sheet_name="Computational")
	no_annotations.to_excel(w, sheet_name="No_annotations")
	not_in_string.to_excel(w, sheet_name="Not_in_string")
	deps_string_ids.to_excel(w, sheet_name="DEPs")





