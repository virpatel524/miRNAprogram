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
mirVIR analyzes enrichment levels for microRNA (miRNA) ages in diseases or any other classes 

Tabs are used to separate the elements of each line.

DEPENDENCIES: NumPy, matplotlib, DendropPy, and SciPy
Pip (https://pypi.python.org/pypi/pip) provides the simplest way to install each dependency

Command line arguments:

-h: display this message

-a: file of miRNAs with ages in the following format:
miRNAname	age
Ages should be quantitative, not expressed in terms of pylogenetic tree distance
To obtain an age file for your set of miRNAs, we recommend you use ProteinHistorian (doi:10.1371/journal.pcbi.1002567)

-f: list of miRNA families and members of the families in the following format for each line:
fam1 mem1
fam1 mem2
fam2 mem1
fam1 mem3

for example

-d: file containing a miRNA with its disease or other biological function
mirna associated_disease
UPREGULATE DOWNREGULATE????

-t: phylogenetic tree file in the Newick format. Can be used to determine age enrichments for certain clades. This currently does nothing.

-y: specifies the format of the ages presented in the age file. If the user gives a parameter that starts with the letter m ('mya', 'million-years-ago',
	'muffin'), then the program interprets the age readings as actual ages rather than distances on a phylogenetic tree

If you wish for the miRNA feature you study to be miRNA families, simply replace the disease file with the family file. If you would like 
the complte age analyses, please repeat using the -f command


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
	
def treelables(treefile):

	fle = open(treefile,"r")
	text = fle.readlines()
	fle.close()
	labeldict = {}

	for line in text:
		if line[0] == "#":
			continue
		pie = string.replace(line, "\n", "")
		pie = pie.split("\t")
		labeldict[float(pie[-1])] = pie[0]
	return labeldict





def makeplot(mirna2age,disfle,treefle=""):
	

	taxdict = treelables(treefle)


	data = mirna2age.values()
	alldata = sorted(data)
	labs = sorted(list(set(data)))

	labels = []

	for i in range(len(labs)):
		labels.append(taxdict[float(labs[i])])

	howmany = {}
	for el in labs:
		howmany[el] = 0

	for el in alldata:
		howmany[el] += 1





	bins = []

	for el in labs:
		bins.append(howmany[el])

	left = range(len(bins))

	plt.bar(left, bins)
	plt.xticks(left,labels,rotation=35)
	plt.title("miRViewer Ages Bincount - New Fam Definition")
	plt.xlabel("Frequency of Ages")





	if not os.path.exists("bincounts"):
		os.makedirs("bincounts")	
	if not os.path.exists("txtfiles"):
		os.makedirs("txtfiles")

	
	os.chdir("bincounts")

	plt.savefig("agesbincount.png", bbox_inches='tight')

	return

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
	if not os.path.exists("enrichcomp"):
		os.makedirs("enrichcomp")
	os.chdir("enrichcomp")


	allages = sorted(mirna2age.values())
	maxage,minage = max(allages), min(allages)

	bgsp = [0] * int(maxage + 1)

	for index, val in enumerate(bgsp):
		bgsp[index] = float(agenumber(index,allages)) / float(len(allages))


	os.chdir("..")

	os.chdir("txtfiles")

	fle = open("disvsbg.txt","w")

	PREAMBLE = "CLASS\tNUM MEM\tAVG AGE\tMED AGE\tPVALUE\n"
	PREAMBLE += "Background\t%i\t%.4f\t%.4f\t\n" %(len(allages),np.mean(allages),np.median(allages))

	fle.write(PREAMBLE)

	os.chdir("..")

	os.chdir("enrichcomp")

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
			ax.set_title('miRNA Age Distribution Comparison \n' + mwu_str)
			ax.set_xticks(ind)

			ymin = ax.viewLim.ymin
			ymax = ax.viewLim.ymax
			y_range = ymax - ymin
			ax.set_ylim(-.05 * y_range, ymax + (.23 * y_range))
			
			ax.legend( (rects1[0], rects2[0]), ('Background (N = ' + str(len(mirna2age)) +")" , 'Class (N = ' +str(len( dis2mirna[dis])  ) + ")") ) 

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

		temp1 = validcheck(temp[1],posmirna);
		if len(temp1) > 0:
			fam2kids.setdefault(temp[0],[]).append(temp1)

		for el in validcheck(temp[1],posmirna):
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


def classagetesting(dis2mirna,kids2age):
	print "CLASS AGE ANALYSIS"
	clsaverage = {}
	clsstd = {}
	clsramavg = {}
	clsramstd = {}
	for cls in dis2mirna:
		if len(dis2mirna[cls]) > 1 :
			clsaverage[cls] = np.mean(dictmap(kids2age, dis2mirna[cls]))
			clsstd[cls] = np.std(dictmap(kids2age, dis2mirna[cls]))

	number = 0
	clsnum = len(clsstd.keys())		
	for cls in dis2mirna:
		if len(dis2mirna[cls]) > 1 :
			clsramavg[cls],clsramstd[cls] = randomages(kids2age, len(dis2mirna[cls]))
			number += 1
			print "CLASS ",number, "out of ", clsnum, "completed"


	avgp = {}
	stdp = {}

	for cls in clsramavg:
		counter = 0
		count2 = 0

		for el in clsramavg[cls]:
			if abs(float(clsaverage[cls]) - float(el))  < .0001:
				counter += 1
		for var in clsramstd[cls]:
			if abs(float(clsstd[cls]) - float(var)) < .0001:
				count2 += 1



		avgp[cls] = float(counter) / float(10000)
		stdp[cls] = float(count2) / float(10000)




	fle = open("clsavage.txt","w")

	fle.write("CLASS"+"\t"+"# of mem"+ "\t"+"avg age"+"\t"+"avg p-value" + "\n")

	for cls in avgp:
		
		fle.write("%s\t%i\t%.2f\t%.4e\n" %(cls, len(dis2mirna[cls]), clsaverage[cls], avgp[cls]))


	fle = open("clsstdage.txt","w")

	fle.write("CLASS"+"\t"+"# of mem"+ "\t"+"std of age"+"\t"+"std p-value" + "\n")

	for cls in stdp:
		fle.write("%s\t%i\t%.2f\t%.4e\n" %(cls, len(dis2mirna[cls]), clsstd[cls], stdp[cls]))

	return 


# Correlation between the age of a disease miRNA and the number of diseases it's associated with.
# More correlations can be added later



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


	fle = open("enrichagetxt.txt","w")

	fle.write("Correlation between the age of a miRNA and the number of classes it is associated with\n")
	fle.write("Coefficient: %.4f\n" % (agedislen))
	fle.write("p-value: %.4e\n" % (agedis_pval))
	fle.close()

	plt.close()
	

	plt.title("miRNA-Class Association Relationships\nCorrelation between the age of a miRNA and the number of classes it is associated with: %.4f\n P-Value: %.4e\n " %(agedislen,agedis_pval))
	plt.xlabel("Age Bins") 
	plt.ylabel('Number of Class Associations')

	



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


def  graphs(mirna2age,mirna2dis,dis2mirna):
	mirnaages = []
	numdis = []
	for mirna in mirna2dis:
		if mirna not in mirna2age:
			continue
		else:
			mirnaages.append(mirna2age[mirna])
			numdis.append(len(mirna2dis[mirna]))



	plt.scatter(mirnaages, numdis)
	plt.title("The Age of a of a microRNA Versus Number of Associated Diseases\n.693 Correlation (p-value: 1.8e-75)")
	plt.xlabel("Age of microRNAs")
	plt.ylabel("Number of Diseases for the microRNA")

	plt.savefig("agevsnumberofdis.png")
	






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
	myaorother = "v"
	print opts

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
	if treefle == "":
		sys.exit("PROGRAM WORKS BEST WITH TREE FILE")


	#			
		
	mirna2age,age2mirna = mirnaagedicts(agefle)
	
	mirnawage = mirna2age.keys()
	makeplot(mirna2age,enrichfle,treefle)
	dis2mirna, mirna2dis = diseasedict(enrichfle,mirna2age)
	
	mirnawdis = dis2mirna.keys()
	posmirna = (mirnawage,mirnawdis)
	
	
	fam2kids = {}
	kids2fam = {}

	if famfle != "":
		fam2kids, kids2fam = famdicts(famfle,posmirna)

	dis2age, age2dis, mirnanodis = disagedict(dis2mirna, mirna2age, mirna2dis)
	



	# classagetesting(dis2mirna, mirna2age)


	print "--------------------------------------------------"
	print "AGE ENRICHMENTS"

	print "NOTE: IMAGES ONLY GENERATED FOR CLASSES WITH 4 OR MORE ASSOCIATED MIRNAS"
	# time.sleep(2)


	# makenrichplots(dis2mirna,mirna2age)

	print "PLOT GENERATION COMPLETED"

	othercorrs(mirna2age,mirna2dis,dis2mirna)

	print "ALL ANALYSES COMPLETED"

	if myaorother[0] == "m":
		milyearanalysis(mirna2age)


	# minianal(mirna2age,dis2mirna)


	graphs(mirna2age,mirna2dis,dis2mirna)





main()