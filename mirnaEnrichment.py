import sys
import getopt
import os
import string
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.stats as scs
import time
import matplotlib.axes as axis
import pylab

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

-t: phylogenetic tree file in the Newick format. Can be used to determine age enrichments for certain clades. This currently does nothing.


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

def select(lst, setting):
	newlst = []

	for el in lst:
		if eval(setting):
			newlst.append(el)

	return el
	



def makeplot(mirna2age,agefle):
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
	plt.xlabel("Ages")
	plt.ylabel("Frequency")
	plt.title("miRNA Age Frequencies")
	if not os.path.exists("results"):
		os.makedirs("results")
	os.chdir("results")

	filename = agefle
	filename  = filename.split("/")

	filename = filename[-1]

	if not os.path.exists(filename):
		os.makedirs(filename)
	os.chdir(filename)

	if not os.path.exists("bincounts"):
		os.makedirs("bincounts")	
	if not os.path.exists("txtfiles"):
		os.makedirs("txtfiles")

	
	os.chdir("bincounts")

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


def agenumber(num,lst):
	counter = 0
	for el in lst:
		if el == num:
			counter += 1
	return counter


def enrichfrac(dis2mirna,mirna2age,dis,endpt):
	diseases = dis2mirna[dis]
	ages = sorted(dictmap(mirna2age, diseases))
	perarray = [0] * int(endpt + 1)
	

	for index,val in enumerate(perarray):
		perarray[index] = float(agenumber(index, ages)) / float(len(ages))

	return perarray


# creds to pyplot documentation

def makenrichplots(dis2mirna,mirna2age):
	os.chdir("..")
	if not os.path.exists("diseasecomp"):
		os.makedirs("diseasecomp")
	os.chdir("diseasecomp")


	allages = sorted(mirna2age.values())
	maxage,minage = max(allages), min(allages)

	bgsp = [0] * int(maxage + 1)

	for index, val in enumerate(bgsp):
		bgsp[index] = float(agenumber(index,allages)) / float(len(allages))


	os.chdir("..")

	os.chdir("txtfiles")

	fle = open("disvsbg.txt","w")

	PREAMBLE = "DISEASE\tAVG AGE\tMED AGE\tPVALUE\n"
	PREAMBLE += "Background\t%.4f\t%.4f\t\n" %(np.mean(allages),np.median(allages))

	fle.write(PREAMBLE)

	os.chdir("..")

	os.chdir("diseasecomp")



### YO VIR SHOULD WE ACCOUNT FOR WHEN THERE ARE LESS THAN 4 MIRNAS?

	num4dis = 0

	for dis in dis2mirna:
		if len(dis2mirna[dis]) > 3:
			num4dis += 1

	bgv = {}

	counter = 0

	
	for index,dis in enumerate(dis2mirna):
	
		if len(dis2mirna[dis]) > 3:

			N = maxage
			
			ind = np.arange(N + 1)  # the x locations for the groups
			width = 0.35       # the width of the bars

			fig = plt.figure()
			ax = fig.add_subplot(111)

			
			rects1 = ax.bar(ind, bgsp, width, color='#990000' )

			dsp = enrichfrac(dis2mirna, mirna2age, dis, maxage)
			rects2 = ax.bar(ind+width, dsp, width, color="black")


			mwu, mwu_pval = scs.mannwhitneyu(allages, dictmap(mirna2age,dis2mirna[dis]))
			mwu_str = 'Mann-Whitney U test: U = %.2g (p = %.3g)' % (mwu, mwu_pval)


			# add some
			ax.set_ylabel('Percentages')
			ax.set_title('miRNA Age Enrichments \n' + mwu_str)
			ax.set_xticks(ind)

			

			
			


			ymin = ax.viewLim.ymin
			ymax = ax.viewLim.ymax
			y_range = ymax - ymin
			ax.set_ylim(-.05 * y_range, ymax + (.23 * y_range))
			
			ax.legend( (rects1[0], rects2[0]), ('Background (N = ' + str(len(mirna2age)) +")" , 'Disease (N = ' +str(len( dis2mirna[dis])  ) + ")") ) 

			plt.savefig(dis.replace(" ","") + ".png")
			
			plt.close()

			var1 = np.mean(dictmap(mirna2age,dis2mirna[dis]))

			var2 = np.median(dictmap(mirna2age,dis2mirna[dis]))

			
			fle.write("%s\t%i\t%.4f\t%.4f\t%.4e\n" % (dis, len(dis2mirna[dis]), var1 , var2, mwu_pval))

			counter += 1
			
			print "enrichment plot " + str(counter) + " of " + str(num4dis)

	fle.close()
	return 








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

def diseasedict(disfle,mirna2age):
	fle = open(disfle,"r")
	text = fle.readlines()
	fle.close()

	dis2mirna = {}
	mirna2dis = {}

	for line in text:

		if line[0] == "#": continue

		temp = deletebn(line)
		temp = splitt(temp)

		if temp[0] not in mirna2age:
			continue

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
	avage = []
	stdage = []
	for i in range(0,10000):
		ages = random.sample(allages, length)


		avage.append(np.mean(ages))
		stdage.append(np.std(ages))

	return avage,stdage


def famagetesting(fam2kids,kids2age):
	print "FAMILY AGE ANALYSIS"
	famaverage = {}
	famstd = {}
	famramavg = {}
	famramstd = {}
	for fam in fam2kids:
		if len(fam2kids[fam]) > 1 :
			famaverage[fam] = np.mean(dictmap(kids2age, fam2kids[fam]))
			famstd[fam] = np.std(dictmap(kids2age, fam2kids[fam]))

	number = 0
	famnum = len(famstd.keys())		
	for fam in fam2kids:
		if len(fam2kids[fam]) > 1 :
			famramavg[fam],famramstd[fam] = randomages(kids2age, len(fam2kids[fam]))
			number += 1
			print "family ",number, "out of ", famnum, "completed"





	

	avgp = {}
	stdp = {}

	for fam in famramavg:
		counter = 0
		count2 = 0

		for el in famramavg[fam]:
			if abs(float(famaverage[fam]) - float(el))  < .0001:
				counter += 1
		for var in famramstd[fam]:
			if abs(float(famstd[fam]) - float(var)) < .0001:
				count2 += 1



		avgp[fam] = float(counter) / float(10000)
		stdp[fam] = float(count2) / float(10000)




	fle = open("famavage.txt","w")

	fle.write("family"+"\t"+"# of mem"+ "\t"+"avg age"+"\t"+"avg p-value" + "\n")

	for fam in avgp:
		fle.write(str(str(fam) + "\t" + str(len(fam2kids[fam])) +"\t"+ str(round(famaverage[fam],4))  + "\t"  + str(round(avgp[fam],4))   + "\n"))

	fle.close()
	fle = open("famstdage.txt","w")

	fle.write("family"+"\t"+"# of mem"+ "\t"+"std of age"+"\t"+"std p-value" + "\n")

	for fam in stdp:
		fle.write(str(str(fam) + "\t" + str(len(fam2kids[fam])) +"\t"+ str(round(famstd[fam],4)) +"\t" + str(round(stdp[fam],4))   + "\n"))

	return 	


# Determine  correlation:

# Correlation between the age of a disease miRNA and the number of diseases it's associated with



def othercorrs(mirna2age,mirna2dis,dis2mirna):
	ageslst1 = []
	numdislst1 = []
	enum = []
	for mirna in mirna2dis:
		ageslst1.append(mirna2age[mirna])
		numdislst1.append(len(mirna2dis[mirna]))
		enum.append((int(mirna2age[mirna]),  len(mirna2dis[mirna]) ))

	agedislen,agedis_pval = scs.spearmanr(ageslst1,numdislst1)

	agebins = []

	for i in range(0,int(max(mirna2age.values()) + 1)):
		agebins.append([])

	numberlst = range(0, int(max(mirna2age.values())+1))

	for age,length in enum:
		agebins[age].append(length)


	fle = open("diseaseagesize.txt","w")

	fle.write("Correlation between the age of a miRNA and the number of diseases it is associated with\n")
	fle.write("Coefficient: %.4f\n" % (agedislen))
	fle.write("p-value: %.4e\n" % (agedis_pval))
	fle.close()

	plt.close()
	

	plt.title("miRNA-Disease Association Relationships")
	plt.xlabel("Age Bins") 
	plt.ylabel('Number of Disease Associations')

	plt.boxplot(agebins)

	pylab.xticks([x+1 for x in numberlst],numberlst)
	

	os.chdir("..")

	if not os.path.exists("diseasecor"):
		os.makedirs("diseasecor")

	os.chdir("diseasecor")

	plt.savefig("mirnadisrel.png")
	plt.close()




	return


def milyearanalysis(mirna2age):
	bincount = list(np.bincount(mirna2age.values()))

	years = sorted( list( set( mirna2age.values())))

	agematrix = [0.0] * (len(years) - 1)

	for index,age in enumerate(agematrix):
		agematrix[index] = float(bincount[index]) / ( float(years[index+1]) - float(years[index]))

	

	os.chdir("..")

	os.chdir("txtfiles")

	fle = open("rateofacq.txt","w")
	fle.write("Average Rates of miRNA Acquisition\n")

	for index,el in enumerate(agematrix):
		fle.write("%.2f new miRNAs per million years between %.2f MYA and %.2f MYA\n" %(el, years[index], years[index+1] )  )

	return 



















def main():

    # take command line arguments and put them in an interpretable format
	try:
		opts, args = getopt.getopt(sys.argv[1:], "a:e:h:d:f:t:y:")
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
	myaorother = ""

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
		if el[0] == "-y":
			myaorother = el[1]


	if agefle == "":
		sys.exit("YOU MUST HAVE AN AGE FILE")
	if enrichfle == "":
		sys.exit("YOU MUST HAVE AN ASSOCIATION FILE")		
		
	mirna2age,age2mirna = mirnaagedicts(agefle)
	
	mirnawage = mirna2age.keys()
	makeplot(mirna2age,agefle)
	dis2mirna, mirna2dis = diseasedict(enrichfle,mirna2age)
	
	mirnawdis = dis2mirna.keys()
	posmirna = (mirnawage,mirnawdis)
	
	
	fam2kids = {}
	kids2fam = {}

	if famfle != "":
		fam2kids, kids2fam = famdicts(famfle,posmirna)

	dis2age, age2dis, mirnanodis = disagedict(dis2mirna, mirna2age, mirna2dis)
	



	if famfle != "":
		# famagetesting(fam2kids, mirna2age)
		print "FAMILY AGE ANALYSIS COMPLETED"


	print "--------------------------------------------------"
	print "DISEASE AGE ENRICHMENT"

	print "NOTE: IMAGES ONLY GENERATED FOR DISEASES WITH MORE THAN 4 ASSOCIATED MIRNAS"
	time.sleep(2)


	# makenrichplots(dis2mirna,mirna2age)

	print "PLOT GENERATION COMPLETED"

	othercorrs(mirna2age,mirna2dis,dis2mirna)

	print "ALL ANALYSES COMPLETED"

	if myaorother[0] == "m":
		milyearanalysis(mirna2age)


	








	










main()