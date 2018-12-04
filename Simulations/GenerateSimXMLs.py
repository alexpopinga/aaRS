import os
import shutil

x = 6
y = 6
z = 6
sequenceData = []
simData = []
trees = []

for x in range(6,101):
	with open("/Users/apop146/Downloads/Simulations/IterativeSubstitutionMatrix_sampleTemplate.xml", "r") as templateXML:
		newXML = open("IterativeSubstitutionMatrix_sim" + str(x) + ".xml", "w")
		newXML.write(templateXML.read())
	x += 1
	newXML.close()
	templateXML.close()

for y in range(6,101):
	with open("ClassI_sim" + str(y) + ".fasta", "r") as fastaFile:
		while True:
			data = fastaFile.readline()
			if ">" not in data:
				sequenceData.append(data)
			if not data:
				break
	#print os.path.basename(fastaFile.name)
	fastaFile.close()
	simData.append(sequenceData)
	sequenceData = []
	y += 1
			
#print len(simData)
#print simData[0]

for z in range(6,101):
	with open("AARS_sim_tree" + str(z) + ".tre", "r") as treeFile:
		tree = treeFile.readline()
		trees.append(tree)
	treeFile.close()
	z += 1
	
print trees
print len(trees)	

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