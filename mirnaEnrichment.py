import sys
import getopt
import os
import string
usage = """
miRich analyzes enrichment levels for microRNA (miRNA) ages in diseases 
or any other biological processes related to miRNAs. 

Tabs are used to seperate the elements of each line

Command line arguments:

-h: print this message

-a: file of miRNAs with ages in the following format:
miRNAname	age
Ages can be in terms of node length in a phylogenetic tree file or it can be age in million years.
To obtain this file, we recommend you use ProteinHistorian (doi:10.1371/journal.pcbi.1002567)

-f: list of miRNA families and members of the family in the following format:
fam_name	mem1	mem2	mem3

-d: file containing a miRNA with its disease or other biological function
mirna associated_disease
UPREGULATE DOWNREGULATE????

-t: phylogenetic tree file in the Newick format. Can be used to determine age enrichments for certain clades

"""

# method that creates a dictionary for a miRNA to its age and all ages to the miRNAs in the respective age categories

def mirnaagedicts(agefileloc):
	fle = open(agefileloc,"r")

	mirna2age = {}
	age2mirna = {}

	text = fle.readlines()
	fle.close()

	for line in text:

		if line[0] == "#": continue
		temp = deletebn(line)
		temp = temp.split("\t")

		mirna2age[temp[0]] = float(temp[1])
		age2mirna.setdefault(float(temp[1]), []).append(temp[0])

	return mirna2age,age2mirna

# method that creates dictionaries to identify families and family members

# helper methods to save time


def deletebn(line):
	return string.replace(line, "\n", "")
def splitt(line):
	return(line.split("\t"))

# method that forms dictionaries identifying which family a miRNA belongs to as well as the children of a miRNA family 

def famdicts(famfle):
	fle = open(famfle,"r")
	text = fle.readlines()
	fle.close()

	fam2kids = {}
	kids2fam = {}

	for line in text:

		if line[0] == "#": continue

		temp = deletebn(line)
		temp = temp.split("\t")

		fam2kids[temp[0]] = temp[1:]

		for el in temp[1:]:
			kids2fam[el] = temp[0]
	return fam2kids, kids2fam

# creates dictionaries identifying the miRNA ages associated with a disease as well as the diseases associated with a miRNA
# also returns a list of all miRNAs that have no disease association according to the data given

def diseasedict(disfle):
	fle = open(disfle,"r")
	text = fle.readlines()
	fle.close()

	dis2mirna = {}
	mirna2dis = {}

	for line in text:

		if line[0] == "#": continue

		temp = deletebn(line)
		temp = splitt(temp)

		dis2mirna.setdefault(temp[1],[]).append(temp[0])
		mirna2dis.setdefault(temp[0],[]).append(temp[1])


	return dis2mirna, mirna2dis

def disagedict(dis2mirna, mirna2age,mirna2dis):
	dis2age = {}
	age2dis = {}
	mirnanodis = []

	for key in dis2mirna:
		temp = dis2mirna[key]
		ageslst = []
		for i in temp:
			if i in mirna2age:
				agelst.append(mirna2age[i])
		dis2age[key] = ageslst

	for mirna in mirna2age:
		if mirna in mirna2dis:
			dis = mirna2dis[mirna]

			for d in dis:
				age2dis.setdefault(mirna2age[mirna2age],[]).appned(d)
		else:
			mirnanodis.append(mirna)

	return dis2age, age2dis, mirnanodis
			










def main():

    # take command line arguments and put them in an interpretable format
	try:
		opts, args = getopt.getopt(sys.argv[1:], "a:e:h:d:f:t:")
	except getopt.GetoptError:
		sys.exit(usage)

	# variable initialization
	mirna2age = {}
	age2mirna = {}
	fam2mirna = {}
	fam2en = {}
	mirna2en = {}

	agefle = ""
	enrichfle = ""
	famfle = ""
	treefle = ""
	for el in opts:
		if el[0] == "-h":
			sys.exit(usage)
		if el[0] == "-a":
			agefle = el[1]
		if el[0] == "-d":
			enrichfle = el[1]
		if el[0] == "-f":
			famfle = el[1]
		if el[0] == "-t":
			treefle = el[1]

	if agefle == "":
		sys.exit("YOU MUST HAVE AN AGE FILE")
	if enrichfle == "":
		sys.exit("YOU MUST HAVE AN ASSOCIATION FILE")		
		
	mirna2age,age2mirna = mirnaagedicts(agefle)

	dis2mirna, mirna2dis = diseasedict(enrichfle)
	
	fam2kids = {}
	kids2fam = {}

	if famfle != "":
		fam2kids, kids2fam = famdicts(famfle)



	dis2age, age2dis, mirnanodis = disagedict(dis2mirna, mirna2age, mirna2dis)x






	


main()