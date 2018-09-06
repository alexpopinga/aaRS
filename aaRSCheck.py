x = 0
with open('AARS.xml') as file:
	line = file.readline()
	y = 1
	#while line:
	for x in range(5, 984):
		sequence = line
		print sequence.count('-')
		line = file.readline()
		x += 1
		y += 1
