f = open("MessyTaxonSet.txt","r")
f1 = open("FixedTaxonSet.txt","w+")

f2 = f.readlines()
for x in f2:
	sep = 'spec="Taxon"/>'
	eliminate = x.split(sep, 1)[0]+'spec="Taxon"/>'
	print eliminate
	f1.write(eliminate)
	
f1.close()