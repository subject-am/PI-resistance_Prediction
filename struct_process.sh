#!/bin/bash
structdir=./1NH0_pdb-files
lovodir=/applic/lovoalign/lovoalign-16.342/bin
reducedir=/applic/reduce/reduce-3.23
mifdir=/applic/IsoMif/IsoMif_150311
isomifdir=/applic/IsoMif/IsoMif_150311
referencedir=.

cd $structdir
for structure in $(ls ./structures/*.pdb | sed s':./structures/::'g | sed s'/.pdb//'g );do
	#Add hydrogen atoms to the protease structure
	$reducedir/reduce -p ./structures/$structure.pdb > ${structure}_H.pdb
	#Align the protease structure with the reference structure
	$lovodir/lovoalign -p1 ${structure}_H.pdb -p2 $referencedir/1NH0_ref_final.pdb -o ${structure}_final.pdb > ${structure}_align.txt
	#Calculate the protease MIFs using the reference cavity
	$mifdir/mif_linux_x86_64 -p ${structure}_final.pdb -g $referencedir/1NH0_ref_final_sph_1.pdb -t ${structure}_mif
	#Comparison between the protease binding site MIF and the reference binding site MIF. The output file contains the Tanimoto coefficient.
	time > ${structure}_isomif_time_grid1.txt $isomifdir/isomif_linux_x86_64 -p1 $referencedir/1NH0_ref_final.mif -p2 ${structure}_mif.mif -s 1 -c 1 -d 1.0 -w -o ./results_control/grid1_d10_
	rm *.pdb
done
