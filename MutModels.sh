#!/bin/bash

while read fasta ; do echo P1P > ${fasta%%.fasta}_pattern ; done < consensus.log
for pattern in $(ls *_pattern) ; do sh pattern-HIVp.sh 1NH0_ref.pdb $pattern ; done

mkdir -p ../1NH0_pdb-files
for folder in $(ls -d 1NH0_ref_*_pattern)
do
	mv ${folder}/$(ls $folder --sort=time | head -n 1) ../1NH0_pdb-files/${folder}.pdb # most recent
	rm -r ${folder}
done
