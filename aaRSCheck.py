x = 0
with open("AARS.xml", 'r') as file:
	for x in range(0,980):
		data = file.readline()
		seq = data[data.find("value")+6:]	
		print (len(seq)-1)
		x += 1
