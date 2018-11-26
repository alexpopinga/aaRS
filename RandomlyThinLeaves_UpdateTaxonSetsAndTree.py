import random

x = 0
y = 0
samples = open("10RandomSamples.xml", "w")
taxonSet = open("10RandomSamples_SampleNames.xml", "w")
traitSet = open("10RandomSamples_TraitSet.xml", "w")
tree = open("10RandomSamples_Tree.xml", "w")

for y in range(0,10):
	with open("IterativeSubMatrix_noKMSKS_20states_SEQUENCES.xml", 'r') as file:
		for x in range(0,961):
			data = file.readline()
			z = random.randint(1,2)
			if (z < 2):
				samples.write(data)
				curName = data[data.find('taxon=')+7:data.find('" totalcount')]
				taxonSet.write('<taxon id="')
				taxonSet.write(curName)
				taxonSet.write('"/>')
				taxonSet.write('\n')
				traitSet.write(curName)
				tree.write(curName)
				tree.write(':0.1, ')
				if 'val_' in curName:
					traitSet.write('=0, ')
				elif 'gly_' in curName:
					traitSet.write('=1, ')
				elif 'ala_' in curName:
					traitSet.write('=2, ')
				elif 'asp_' in curName:
					traitSet.write('=3, ')
				elif 'pro_' in curName:
					traitSet.write('=4, ')
				elif 'glu_' in curName:
					traitSet.write('=5, ')
				elif 'ser_' in curName:
					traitSet.write('=6, ')
				elif 'leu_' in curName:
					traitSet.write('=7, ')
				elif 'arg_' in curName:
					traitSet.write('=8, ')
				elif 'thr_' in curName:
					traitSet.write('=9, ')
				elif 'ile_' in curName:
					traitSet.write('=10, ')
				elif 'gln_' in curName:
					traitSet.write('=11, ')
				elif 'asn_' in curName:
					traitSet.write('=12, ')
				elif 'cys_' in curName:
					traitSet.write('=13, ')
				elif 'tyr_' in curName:
					traitSet.write('=14, ')
				elif 'his_' in curName:
					traitSet.write('=15, ')
				elif 'met_' in curName:
					traitSet.write('=16, ')
				elif 'trp_' in curName:
					traitSet.write('=17, ')
				elif 'lys_' in curName:
					traitSet.write('=18, ')
				elif 'phe_' in curName:
					traitSet.write('=19, ')
				else:
					print "we're broken!"
				traitSet.write('\n')
	y += 1
	if (y < 10):
		samples.write('\n')
		taxonSet.write('\n')
		traitSet.write('\n')
		tree.write('\n')
	file.close()
	x = 1
				
samples.close()
taxonSet.close()
traitSet.close()
tree.close()