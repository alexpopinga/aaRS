import numpy as np
import random

def shuffle(seq, i, j):
	seq = list(seq)
	seq[i], seq[j] = seq[j], seq[i]
	return ''.join(seq)

############ for testing ############
exampleAlignment = ["ABCDEFG","ABCDEFG","ABCDEFG"]
exampleIndices = [0,1,2,3,4,5,6]
shuffledIndices = random.sample(range(0,7),7)	
print shuffledIndices

#sequence
a = 0
#site
b = 0
for a in range(0,3):
	for b in range(0, len(exampleIndices)):
		exampleAlignment[a] = shuffle(exampleAlignment[a], exampleIndices[b], shuffledIndices[b])
		b += 1
	a += 1
print exampleAlignment
############ end testing ############

x = 0
sequences = []
sequenceAlignment = []
newAlignment = []

with open('ClassII.xml', 'r') as file:
	for x in range(0,478):
		for line in file.readlines():
			sequences.append(line[line.find('totalcount=\"20\" value=\"')+23:line.find('\"/>')])
			x += 1
print sequences

indexList = []

y = 0
for y in range(0,len(sequences[0])):
	indexList.append(y)
	y += 1
	
orderedIndices = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104]
random.shuffle(indexList)
print indexList

seqNum = 0
siteNum = 0
for seqNum in range(0, len(sequences)):
	for siteNum in range(0, len(orderedIndices)):
		sequences[seqNum] = shuffle(sequences[seqNum], orderedIndices[siteNum], indexList[siteNum])
		#print orderedIndices[siteNum]
		#print indexList[siteNum]
		siteNum += 1
	seqNum += 1
print sequences
