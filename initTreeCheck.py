x = 0
seqCount = 0
taxonCount = 0
taxonNames = []
treeNames = []
def findnth(string, substring, n):
    parts = string.split(substring, n + 1)
    if len(parts) <= n + 1:
        return -1
    return len(string) - len(parts[-1]) - len(substring)
with open('AARS.xml', 'r') as file:
	for line in file.readlines():
		for word in line.split():
			if word == '<sequence':
				seqCount = seqCount + 1
				taxonNames.append(line[line.find('taxon=')+7:line.find('" totalcount')])
			elif word == 'newick=':
				treeNames.append(line[line.find('newick=')+15:line.find(':0.1,')])
				x += 1
				currentTree = line[line.find('newick=')+15:line.find(':0.1):100):100):300);')]
				#print currentTree
				for x in range(1,seqCount):
					print findnth(currentTree, treeNames[x-1], x)
					treeNames.append(line[line.find(':0.1,')+1:line.find(':0.1,')])
					x += 1
		#else:
		#	currentTree = line[line.find(':0.1,')+1:]
			#print currentTree
			#treeNames.append(line[line.find(':0.1,')+1:line.find('')])
		#	taxonCount += 1
		#x += 1
#print treeNames
#print currentTree