import sys
import getopt
import os
import string
import matplotlib.pyplot as plt
import numpy as np
import scipy as scs
import scipy.stats as sstats
import random

usage = """
miRich analyzes enrichment levels for microRNA (miRNA) ages in diseases 
or any other biological processes related to miRNAs. 

Tabs are used to seperate the elements of each line

Command line arguments:

-h: display this message

-a: file of miRNAs with ages in the following format:
miRNAname	age
Ages should be quantitative, not expressed in terms of pylogenetic tree distance
To obtain an age file for your set of miRNAs, we recommend you use ProteinHistorian (doi:10.1371/journal.pcbi.1002567)

-f: list of miRNA families and members of the families in the following format for each line:
fam	mem1	mem2	mem3
fam2	mem1	mem2	mem3

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
	if not os.path.exists("bincounts"):
		os.makedirs("bincounts")	
	if not os.path.exists("txtfiles"):
		os.makedirs("txtfiles")

	os.chdir("..")
	os.chdir("results/bincounts")

	plt.savefig("agesbincount.png")

	os.chdir("..")
	os.chdir("txtfiles")

	agesbin = {}

	for el in mirna2age.values():
		if el not in agesbin:
			agesbin[el] = 1
		else:
			agesbin[el] += 1

	fle = open("agecounts.txt","w")

	for key in agesbin:
		fle.write(str(key) + "\t" + str(agesbin[key]) + "\n")
	
	plt.close()

def validcheck(el,posmirna):
	returner = []
	for i in el:
		if i in posmirna[0] or i in posmirna[1]:
			returner.append(i)
	return returner

# method that forms dictionaries identifying which family a miRNA belongs to as well as the children of a miRNA family
# we have to account for miRNAs in the family file for which there are no ages. As such, in fam2kids, we do not include  

def famdicts(famfle,posmirna):
	fle = open(famfle,"r")
	text = fle.readlines()
	fle.close()

	fam2kids = {}
	kids2fam = {}

	for line in text:

		if line[0] == "#": continue

		temp = deletebn(line)
		temp = temp.split("\t")

		temp1 = validcheck(temp[1:],posmirna);
		if len(temp1) > 0:
			fam2kids[temp[0]] = temp1

		for el in validcheck(temp[1:],posmirna):
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
			
def dictmap(dic, lst):
	newlst = []
	for el in lst:
		newlst.append(dic[el])
	return newlst

def randomages(kids2age,length):
	allages = kids2age.values()
	agedicts = []
	for i in range(0,10000):
		ages = random.sample(allages, length)

		agedicts.append(np.mean(ages))
	return agedicts


def famagetesting(fam2kids,kids2age):
	famaverage = {}
	fammedian = {}
	famstd = {}
	fam2random = {}
	for fam in fam2kids:
		if fam2kids[fam] != []:
			famaverage[fam] = np.mean(dictmap(kids2age, fam2kids[fam]))
			fammedian[fam] = np.median(dictmap(kids2age, fam2kids[fam]))

	number = 0
	famnum = len(fam2kids.keys())		
	for fam in fam2kids:
		fam2random[fam] = randomages(kids2age, len(fam2kids[fam]))
		number += 1
		print "family completed:", number, "out of ", famnum

	



	


	







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
	
	mirnawage = mirna2age.keys()
	makeplot(mirna2age)
	dis2mirna, mirna2dis = diseasedict(enrichfle)
	
	mirnawdis = dis2mirna.keys()
	posmirna = (mirnawage,mirnawdis)
	
	
	fam2kids = {}
	kids2fam = {}

	if famfle != "":
		fam2kids, kids2fam = famdicts(famfle,posmirna)

	dis2age, age2dis, mirnanodis = disagedict(dis2mirna, mirna2age, mirna2dis)
	famagetesting(fam2kids, mirna2age)

	








	










main()