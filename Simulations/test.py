x = 0

list1 = ["blah", "bloop", "blorg"]
list2 = ["arf", "meow", "eep"]

listOflists = []
listOflists.append(list1)
listOflists.append(list2)

for x in range(0,3):
	with open("testOutput.xml", "r+") as outfile:
		content = outfile.read()
		outfile.close()
	with open("testOutput.xml", "r+") as outfile:
		outfile.write(content.replace("$BLAH", str(listOflists[0][x]), 1))
		content = outfile.read()
		outfile.close()
		#print content
		#print x