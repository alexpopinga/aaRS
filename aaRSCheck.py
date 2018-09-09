x = 1
with open("AARS.xml", 'r') as file:
	for x in range(1,980):
		data = file.readline()
		seq = data[data.find("value")+6:]
		seqLength = len(seq)-1
		if (seqLength != 116) & (x > 5):
			print("On line ")
			print x
			print(" the sequence length is ")
			print seqLength
		x += 1
