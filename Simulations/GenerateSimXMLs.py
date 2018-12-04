import os
import shutil

x1 = 6
x2 = 6
y = 6
z = 6
sequenceData1 = []
simData1 = []
sequenceData2 = []
simData2 = []
trees = []

for x1 in range(6,101):
	with open("ClassI_sim" + str(x1) + ".fasta", "r") as fastaFile:
		while True:
			data = fastaFile.readline()
			if ">" not in data:
				sequenceData1.append(data)
			if not data:
				break
	#print os.path.basename(fastaFile.name)
	fastaFile.close()
	simData1.append(sequenceData1)
	sequenceData1 = []
	x1 += 1
			
#print len(simData1)
#print simData1[0]

for x2 in range(6,101):
	with open("ClassII_sim" + str(x2) + ".fasta", "r") as fastaFile:
		while True:
			data = fastaFile.readline()
			if ">" not in data:
				sequenceData2.append(data)
			if not data:
				break
	#print os.path.basename(fastaFile.name)
	fastaFile.close()
	simData2.append(sequenceData2)
	sequenceData2 = []
	x1 += 1
	
#print len(simData2)
#print simData2[0]

for y in range(6,101):
	with open("AARS_sim_tree" + str(y) + ".tre", "r") as treeFile:
		tree = treeFile.readline()
		trees.append(tree)
	treeFile.close()
	y += 1
	
#print trees
#print len(trees)	

for z in range(6,101):
	with open("/Users/apop146/Downloads/Simulations/IterativeSubstitutionMatrix_sampleTemplate.xml", "r") as templateXML:
		newXML = open("IterativeSubstitutionMatrix_sim" + str(z) + ".xml", "w")
		newXML.write(templateXML.read())
	z += 1
	newXML.close()
	templateXML.close()
	
a1 = 6
a2 = 0
b = 0

#print len(trees)

for a1 in range(6,101):
	with open("IterativeSubstitutionMatrix_sim" + str(a1) + ".xml", "r") as f1:
		newFile = f1.read()
		f1.close()
		with open("IterativeSubstitutionMatrix_sim" + str(a1) + ".xml", "w") as f2:
			f2.write(newFile.replace("$NEWICKTREE", str(trees[a2])))
			f2.close()
	a1 += 1
	a2 += 1

#print simData1[0][99]	

with open("IterativeSubstitutionMatrix_sim6.xml", "r") as infile:
	content = infile.read()
	infile.close()
	for b in range(0,99):
		with open("IterativeSubstitutionMatrix_sim6.xml", "r+") as outfile:
			outfile.write(content.replace("$SEQUENCEDATA", str(simData1[0][b]), 1))
			#print simData1[0][b]
			content = outfile.read()
			outfile.close()
		b += 1
	
#import random

#x = 0
#y = 0
#fullTaxaList = []
#thinnedTaxaList = []
#newFile = open("ClassICtermClassIINtermSequences.xml", "w")

#with open("thinnedSequences.xml",'r') as thinnedFile:
#	for y in range(0,479):
#		data = thinnedFile.readline()
#		curName = data[data.find('taxon=')+7:data.find('" totalcount')]
#		thinnedTaxaList.append(curName)
#	y += 1
		
#with open("fullSequences.xml", 'r') as fullFile:
#	for x in range(0,961):
#		data = fullFile.readline()
#		curName = data[data.find('taxon=')+7:data.find('" totalcount')]
#		fullTaxaList.append(curName)
#		if fullTaxaList[x] in thinnedTaxaList:
#			newFile.write(data)
#	x += 1		

#fullFile.close()
#thinnedFile.close()				
#newFile.close()