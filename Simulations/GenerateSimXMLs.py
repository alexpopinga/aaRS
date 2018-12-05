import os

x1 = 6
x2 = 6
y = 6
z = 6
sequenceData1 = []
simData1 = []
sequenceData2 = []
simData2 = []
trees = []

########################################
# loop through simulated ClassI data
########################################
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
			
#print len(simData1)
#print simData1[0]

########################################
# loop through simulated ClassII data
########################################
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
	
#print len(simData2)
#print simData2[0]

########################################
# get simulated trees (ClassI and II)
########################################
for y in range(6,101):
	with open("AARS_sim_tree" + str(y) + ".tre", "r") as treeFile:
		tree = treeFile.readline()
		trees.append(tree)
	treeFile.close()
	
#print trees
#print len(trees)	

########################################
# generate 100 new XMLs, 1 for each sim.
########################################
for z in range(6,101):
	with open("/Users/apop146/Downloads/Simulations/IterativeSubstitutionMatrix_sampleTemplate.xml", "r") as templateXML:
		newXML = open("IterativeSubstitutionMatrix_sim" + str(z) + ".xml", "w")
		newXML.write(templateXML.read())
	newXML.close()
	templateXML.close()

########################################
# use simulated trees as init. trees
########################################
a = 6

#print len(trees)

for a in range(6,101):
	with open("IterativeSubstitutionMatrix_sim" + str(a) + ".xml", "r") as f1:
		newFile = f1.read()
		f1.close()
		with open("IterativeSubstitutionMatrix_sim" + str(a) + ".xml", "w") as f2:
			f2.write(newFile.replace("$NEWICKTREE", str(trees[a-6])))
			f2.close()
	
########################################
# place simulated data into XMLs
########################################
b1 = 0
b2 = 0
c = 6

#print simData1[0][99]	

for c in range(6,101):
	for b1 in range(0,100):
		with open("IterativeSubstitutionMatrix_sim" + str(c) + ".xml", "r+") as outfile:
			content = outfile.read()
			outfile.close()
		with open("IterativeSubstitutionMatrix_sim" + str(c) + ".xml", "r+") as outfile:
			outfile.write(content.replace("$SEQUENCEDATA", str(simData1[c-6][b1]), 1))
			#print simData1[0][b1]
			content = outfile.read()
			outfile.close()
	for b2 in range(0,100):
		with open("IterativeSubstitutionMatrix_sim" + str(c) + ".xml", "r+") as outfile:
			content = outfile.read()
			outfile.close()
		with open("IterativeSubstitutionMatrix_sim" + str(c) + ".xml", "r+") as outfile:
			outfile.write(content.replace("$SEQUENCEDATA", str(simData2[c-6][b2]), 1))
			#print simData1[0][b2]
			content = outfile.read()
			outfile.close()