import random

x = 0
y = 0
fullTaxaList = []
thinnedTaxaList = []
newFile = open("ClassICtermClassIINtermSequences.xml", "w")

with open("3sequences.xml",'r') as thinnedFile:
	for y in range(0,3):
		data = thinnedFile.readline()
		curName = data[data.find('taxon=')+7:data.find('" totalcount')]
		thinnedTaxaList.append(curName)
	y += 1
		
with open("5sequences.xml", 'r') as fullFile:
	for x in range(0,5):
		data = fullFile.readline()
		curName = data[data.find('taxon=')+7:data.find('" totalcount')]
		fullTaxaList.append(curName)
		if fullTaxaList[x] in thinnedTaxaList:
			newFile.write(data)
	x += 1		

fullFile.close()
thinnedFile.close()				
newFile.close()