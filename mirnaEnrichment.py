import sys
import getopt
import os
import string
import matplotlib.pyplot as plt
import numpy as np
import scipy as scs

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

def chunks(l, n):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def deletebn(line):
	return string.replace(line, "\n", "")
def splitt(line):
	return(line.split("\t"))
def flatten(l):
	newlst = []
	print l

	for el in l:
		for each in el:
			newlst.append(each)

	return newlst


def boxplots(dis2age,keys):
	lst = []
	for element in keys:
		lst.append(dis2age[element])
	plt.boxplot(lst,keys)
	



def makeplot(mirna2age):
	binlst = map(int, mirna2age.values())
	binlst.sort()
	ranger = range(0,binlst[-1]+1)

	freq = np.bincount(np.array(binlst),weights=None)

	pos = np.arange(len(ranger))
	width = 1.0
	ax = plt.axes()
	ax.set_xticks(pos + (width / 2))
	ax.set_xticklabels(ranger)
	plt.bar(pos,freq,width)
	plt.xlim(pos.min(),pos.max()+width)
	if not os.path.exists("results"):
		os.makedirs("results")
	os.chdir("results")
	plt.savefig("agesbincount.png")

	plt.close()


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

# creates dictionaries mapping diseases to associated miRNAs and vice versa

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

	# we exclude the case where an age is not present in your database. This could be simply due to a lack of data and thus it isn't very helpful

	for key in dis2mirna:
		temp = dis2mirna[key]
		ageslst = []
		for i in temp:
			if i in mirna2age:
				ageslst.append(mirna2age[i])
		dis2age[key] = ageslst

	# here we  account for no disease assocation.	

	for mirna in mirna2age:
		if mirna in mirna2dis:
			dis = mirna2dis[mirna]

			for d in dis:
				age2dis.setdefault(mirna2age[mirna],[]).append(d)
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
	makeplot(mirna2age)
	dis2mirna, mirna2dis = diseasedict(enrichfle)
	
	
	fam2kids = {}
	kids2fam = {}

	if famfle != "":
		fam2kids, kids2fam = famdicts(famfle)

	dis2age, age2dis, mirnanodis = disagedict(dis2mirna, mirna2age, mirna2dis)

	








	











	


main()