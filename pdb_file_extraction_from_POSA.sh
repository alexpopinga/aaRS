#!/bin/bash
i=0
cat Leu_Euk_POSA.pdb > curFile.txt

#Count the number of structures (by "TER") in the POSA alignment
cat curFile.txt | grep -o TER | wc -w > numFiles.txt

read numFiles < numFiles.txt 
while [ $i -lt $numFiles ]; do
	cat Leu_Euk_POSA.pdb | awk '/TER/{j++}j=='$i'' > Leu_Euk_POSA_$i.pdb; echo "TER END" >> Leu_Euk_POSA_$i.pdb
	let i=i+1
done
