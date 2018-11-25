import random

x = 1
y = 1
samples = open("10RandomSamples.xml", "w")

for y in range(1,10):
	with open("IterativeSubMatrix_noKMSKS_20states_SEQUENCES.xml", 'r') as file:
		for x in range(1,961):
			data = file.readline()
			z = random.randint(1,2)
			if (z < 2):
				samples.write(data)
	y += 1
	if (y < 10):
		samples.write("\n")
	file.close()
	x = 1
				
samples.close()