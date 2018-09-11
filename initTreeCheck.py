x = 0
seqCount = 0
taxonCount = 0
taxonNames = []
treeNames = []
curTreeIndex = 0

def findnth(string, substring, n):
    parts = string.split(substring, n+1)
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
				tree = line[line.find('newick=')+15:line.find('):100):100):300);')]
				# kludge-y replace method for parsing
				tree = tree.replace(':0.1, ','|').replace(':0.1):100, (','|').replace(':0.1):100):100, ((','|').replace(':0.1):100):100):100, (((','|').replace(':0.1):100):100):100):100, ((((','|').replace(':0.1):100):100):100):100):100, (((','|').replace(':0.1','|')
				for x in range(0,seqCount):
					newTreeIndex = findnth(tree, '|', x)
					treeNames.append(tree[curTreeIndex:newTreeIndex])
					curTreeIndex = newTreeIndex+1
					# why isn't there "seqCount" number of taxa in tree?
					if treeNames[x] == 'phe_euk_B_hominis':
						if len(taxonNames) != len(treeNames):
							print("These taxa are missing in the initial tree: ")
							print set(taxonNames) - set(treeNames)
							print("These sequence taxa are missing: ")
							print set(treeNames) - set(taxonNames)
						else:
							print treeNames
					x += 1