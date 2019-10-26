#!/bin/bash

# Compare with susceptible
for file in $(ls *.fasta)
do
	for position in {0..98}
	do
		seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
		sequence=$(sed '1d' $file)
		ressus=${seqcons:$position:1}
		resseq=${sequence:$position:1}
		respos=$(expr $position + 1)
		if [[ $ressus != $resseq ]]
		then
			echo $ressus$respos$resseq >> ${file%%.fasta}_pattern
		fi
	done
done
