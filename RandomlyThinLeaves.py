import random

x = 1
y = 1
text_file = open("output.txt", "w")

with open("IterativeSubMatrix_noKMSKS_20states_SEQUENCES.xml", 'r') as file:
	for x in range(1,980):
		data = file.readline()
		z = random.randint(1,2)
		if (z < 2):
			text_file.write(data)
				
text_file.close()