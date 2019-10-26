#!/bin/bash

filename=$1
dataset=${filename%.*}

# Separate sequences w/o info on PI
head -n 1 $dataset > ${dataset}-NA
head -n 1 $dataset > ${dataset}-PI
sed '1d' $dataset | while read line
do
	checkline=$(grep 'NA' <<< $line)
	if [[ -n $checkline ]]
	then
		echo $line >> ${dataset}-NA
	else
		echo $line >> ${dataset}-PI
	fi
done
nlinesNA=$(wc -l ${dataset}-NA | awk '{print $1}')
nlinesPI=$(wc -l ${dataset}-PI | awk '{print $1}')
if [[ $nlinesNA -eq 1 ]]; then rm ${dataset}-NA ; 
elif [[ $nlinesPI -eq 1 ]]; then rm ${dataset}-PI ; echo "Review script parameters for susceptibility." ; exit 1 ; fi
if [[ -f ${dataset}-NA ]] ; then sed -i '/^$/d' ${dataset}-NA ; echo "The sequences which don't have information on PIs are saved in "${dataset}-NA ; fi
sed -i '/^$/d' ${dataset}-PI ; echo "The sequences with information on resistance to PIs are saved in "${dataset}-PI

#### Susceptible to all PI ####
head -n 1 $dataset > ${dataset}-susceptible
head -n 1 $dataset > ${dataset}-res20
sed '1d' ${dataset}-PI | while read sequence
do
	DRV=$(awk '{print $9}' <<< $sequence) ; sDRV=${DRV%%.*} # don't use printf %1.0f because of approximations

	FPV=$(awk '{print $2}' <<< $sequence) ; sFPV=${FPV%%.*}
	ATV=$(awk '{print $3}' <<< $sequence) ; sATV=${ATV%%.*}
	IDV=$(awk '{print $4}' <<< $sequence) ; sIDV=${IDV%%.*}
	LPV=$(awk '{print $5}' <<< $sequence) ; sLPV=${LPV%%.*}
	NFV=$(awk '{print $6}' <<< $sequence) ; sNFV=${NFV%%.*}
	SQV=$(awk '{print $7}' <<< $sequence) ; sSQV=${SQV%%.*}
	TPV=$(awk '{print $8}' <<< $sequence) ; sTPV=${TPV%%.*}

	if [[ $sFPV -lt 3 && $sATV -lt 3 && $sIDV -lt 3 && $sLPV -lt 3 && $sNFV -lt 3 && $sSQV -lt 3 && $sTPV -lt 3 && $sDRV -lt 3 ]]
	then
		echo $sequence >> ${dataset}-susceptible
	elif [[ $sFPV -gt 20 && $sATV -gt 20 && $sIDV -gt 20 && $sLPV -gt 20 && $sNFV -gt 20 && $sSQV -gt 20 && $sTPV -gt 20 && $sDRV -gt 20 ]]
	then
		echo $sequence >> ${dataset}-res20
		echo $sequence >> ${dataset}-res15
		echo $sequence >> ${dataset}-res10
	elif [[ $sFPV -gt 15 && $sATV -gt 15 && $sIDV -gt 15 && $sLPV -gt 15 && $sNFV -gt 15 && $sSQV -gt 15 && $sTPV -gt 15 && $sDRV -gt 15 ]]
	then
		echo $sequence >> ${dataset}-res15
		echo $sequence >> ${dataset}-res10
	elif [[ $sFPV -gt 10 && $sATV -gt 10 && $sIDV -gt 10 && $sLPV -gt 10 && $sNFV -gt 10 && $sSQV -gt 10 && $sTPV -gt 10 && $sDRV -gt 10 ]]
	then
		echo $sequence >> ${dataset}-res10
	fi
done
nlinesSUS=$(wc -l ${dataset}-susceptible | awk '{print $1}')
if [[ $nlinesSUS -eq 1 ]]; then rm ${dataset}-susceptible ; echo "No sequence is susceptible to all PI. Review script parameters for susceptibility." ; exit 0 ; fi
sed -i '/^$/d' ${dataset}-susceptible ; echo "All sequences susceptible to DRV, FPV, ATV, IDV, LPV, NFV, SQV, and TPV are saved in "${dataset}-susceptible
nlinesRES20=$(wc -l ${dataset}-res20 | awk '{print $1}')
if [[ $nlinesRES20 -eq 1 ]]; then rm ${dataset}-res20 ; echo "No sequence is highly resistant to all PI. Review script parameters for susceptibility." ; exit 0 ; fi
sed -i '/^$/d' ${dataset}-res20 ; echo "All sequences highly resistant to DRV, FPV, ATV, IDV, LPV, NFV, SQV, and TPV are saved in "${dataset}-res20
nlinesRES15=$(wc -l ${dataset}-res15 | awk '{print $1}')
if [[ $nlinesRES15 -eq 1 ]]; then rm ${dataset}-res15 ; echo "No sequence is highly resistant to all PI. Review script parameters for susceptibility." ; exit 0 ; fi
sed -i '/^$/d' ${dataset}-res15 ; echo "All sequences highly resistant to DRV, FPV, ATV, IDV, LPV, NFV, SQV, and TPV are saved in "${dataset}-res15


#### Getting the sequences ####
sed '1d' ${dataset}-susceptible | while read sequence
do
	seqID=$(awk '{print $1}' <<< $sequence)
	echo ">"$seqID > ${seqID}.fasta
	checkseqID=$(ls | grep $seqID | wc -l)
	
	seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
	# consensus B subtype (https://hivdb.stanford.edu/pages/documentPage/consensus_amino_acid_sequences.html)
	
	for columnno in {10..108}
	do
		aminoacid=$(awk '{print $'$columnno'}' <<< $sequence)
		position=$(expr $columnno - 10)
		
		rescon=${seqcons:$position:1}
		respos=$(expr $position + 1)
		
		checkaa=$(wc -m <<< $aminoacid) # in dataset, if there's a mixture the $(wc -m) is ge 3, else is 2
		if [[ -n $checkseqID && $checkaa -eq 2 ]]
		then
			for file in $(ls | grep $seqID)
			do
				if [[ $aminoacid == '-' ]] ; then printf "%s" $rescon >> $file
				elif [[ $aminoacid == '.' || $aminoacid == '~' ]] ; then
					echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}d.fasta
                    printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}d.fasta
				else printf "%s" $aminoacid >> $file ; fi
			done
		elif [[ -n $checkseqID && $checkaa -ge 3 ]]
		then
			for file in $(ls | grep $seqID | grep -v STOP)
			do
				aa=$(expr $checkaa - 2) # discount for (newline + count from zero)
				for resmix in $(seq 0 $aa)
				do
					resaa=${aminoacid:$resmix:1} # keep mixtures with consensus aminoacid, because of further mutations
					if [[ -n $(grep '*' <<< $resaa) ]] ; then
						echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}_STOP.fasta
						printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}_STOP.fasta
					elif [[ -n $(grep '#' <<< $resaa) ]] ; then
	                    echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
						printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
                        printf "%s" $aminoacid >> ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
					else
						echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
                        printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
                        printf "%s" $resaa >> ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
					fi
				done
				rm $file
			done
		fi
	done
done
echo "All susceptible sequences separated in FASTA files."

# Check for X residue
sed '1d' ${dataset}-susceptible | while read sequence
do
	seqID=$(awk '{print $1}' <<< $sequence)
	for file in $(ls | grep $seqID)
	do
		checkX=$(sed '1d' $file | grep X)
		if [[ -n $checkX ]] ; then echo $seqID >> ${dataset}-susceptible_Xmut.tmp ; fi
	done
done
if [[ -f ${dataset}-susceptible_Xmut.tmp ]]
then
	cat PI_dataset-susceptible_Xmut.tmp | sort | uniq > PI_dataset-susceptible_Xmut.log
	rm PI_dataset-susceptible_Xmut.tmp
	echo "The following sequences have at least one mutation for any possible amino acid:"
	echo
	cat ${dataset}-susceptible_Xmut.log
	echo
	echo "The IDs are saved in the "${dataset}-susceptible_Xmut.log" file."

# Separate various X mutations
	while read xmut
	do
		for file in $(ls | grep $xmut)
		do
			seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
			seqX=$(sed '1d' $file | grep X)
			seqNX=$(sed '1d' $file | grep -o X | wc -l)
			if [[ $seqNX -eq 1 ]]
			then
				for position in {0..98}
				do
					resseq=${seqX:$position:1}
					rescon=${seqcons:$position:1}
					respos=$(expr $position + 1)
				
					## Amino acids
					# "A" = "ALA"
					# "R" = "ARG"
					# "D" = "ASP"
					# "N" = "ASN"
					# "C" = "CYS"
					# "E" = "GLU"
					# "Q" = "GLN"
					# "G" = "GLY"
					# "H" = "HIS"
					# "I" = "ILE"
					# "L" = "LEU"
					# "K" = "LYS"
					# "M" = "MET"
					# "F" = "PHE"
					# "P" = "PRO"
					# "S" = "SER"
					# "T" = "THR"
					# "W" = "TRP"
					# "Y" = "TYR"
					# "V" = "VAL"
			
					if [[ $resseq == "X" ]]
					then
						for aa in A R D N C E Q G H I L K M F P S T W Y V
						do
							echo ">"$xmut > ${file%%.fasta}_${rescon}${respos}${aa}.fasta
							sed -r 's/(.{'$position'})(.)/\1'$aa'/' <<< $seqX >> ${file%%.fasta}_${rescon}${respos}${aa}.fasta
						done
					fi
				done
				rm $file
				echo "All sequences from X mutations available."
			fi
		done
	done < ${dataset}-susceptible_Xmut.log
fi

# Move FASTA to folder
rundate=$(date "+%Y%m%d")
mkdir -p ${rundate}_${dataset}_susceptible
mv *.fasta ${rundate}_${dataset}_susceptible

sed '1d' ${dataset}-res20 | while read sequence
do
	seqID=$(awk '{print $1}' <<< $sequence)
	echo ">"$seqID > ${seqID}.fasta
	checkseqID=$(ls | grep $seqID | wc -l)
	
	seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
	# consensus B subtype (https://hivdb.stanford.edu/pages/documentPage/consensus_amino_acid_sequences.html)
	
	for columnno in {10..108}
	do
		aminoacid=$(awk '{print $'$columnno'}' <<< $sequence)
		position=$(expr $columnno - 10)
		
		rescon=${seqcons:$position:1}
		respos=$(expr $position + 1)
		
		checkaa=$(wc -m <<< $aminoacid) # in dataset, if there's a mixture the $(wc -m) is ge 3, else is 2
		if [[ -n $checkseqID && $checkaa -eq 2 ]]
		then
			for file in $(ls | grep $seqID)
			do
				if [[ $aminoacid == '-' ]] ; then printf "%s" $rescon >> $file
				elif [[ $aminoacid == '.' || $aminoacid == '~' ]] ; then
					echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}d.fasta
					printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}d.fasta
				else printf "%s" $aminoacid >> $file ; fi
			done
		elif [[ -n $checkseqID && $checkaa -ge 3 ]]
		then
			for file in $(ls | grep $seqID)
			do
				aa=$(expr $checkaa - 2) # discount for (newline + count from zero)
				for resmix in $(seq 0 $aa)
				do
					resaa=${aminoacid:$resmix:1} # keep mixtures with consensus aminoacid, because of further mutations
					if [[ -n $(grep '*' <<< $resaa) ]] ; then
                                                echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}_STOP.fasta
	                                        printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}_STOP.fasta
					elif [[ -n $(grep '#' <<< $resaa) ]] ; then 
	        				echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
		   				printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
			                        printf "%s" $aminoacid >> ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
                                        else
						echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
						printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
						printf "%s" $resaa >> ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
					fi
				done
				rm $file
			done
		fi
	done
done
echo "All highly resistant sequences separated in FASTA files."

# Check for X residue
sed '1d' ${dataset}-res20 | while read sequence
do
	seqID=$(awk '{print $1}' <<< $sequence)
	for file in $(ls | grep $seqID)
	do
		checkX=$(sed '1d' $file | grep X)
		if [[ -n $checkX ]] ; then echo $seqID >> ${dataset}-res20_Xmut.tmp ; fi
	done
done
if [[ -f ${dataset}-res20_Xmut.tmp ]]
then
	cat PI_dataset-res20_Xmut.tmp | sort | uniq > PI_dataset-res20_Xmut.log
	rm PI_dataset-res20_Xmut.tmp
	echo "The following sequences have at least one mutation for any possible amino acid:"
	echo
	cat ${dataset}-res20_Xmut.log
	echo
	echo "The IDs are saved in the "${dataset}-res20_Xmut.log" file."

# Separate various X mutations
	while read xmut
	do
		for file in $(ls | grep $xmut)
		do
			seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
			seqX=$(sed '1d' $file | grep X)
			seqNX=$(sed '1d' $file | grep -o X | wc -l)
			if [[ $seqNX -eq 1 ]]
			then
				for position in {0..98}
				do
					resseq=${seqX:$position:1}
					rescon=${seqcons:$position:1}
					respos=$(expr $position + 1)

					## Amino acids
					# "A" = "ALA"
					# "R" = "ARG"
					# "D" = "ASP"
					# "N" = "ASN"
					# "C" = "CYS"
					# "E" = "GLU"
					# "Q" = "GLN"
					# "G" = "GLY"
					# "H" = "HIS"
					# "I" = "ILE"
					# "L" = "LEU"
					# "K" = "LYS"
					# "M" = "MET"
					# "F" = "PHE"
					# "P" = "PRO"
					# "S" = "SER"
					# "T" = "THR"
					# "W" = "TRP"
					# "Y" = "TYR"
					# "V" = "VAL"

					if [[ $resseq == "X" ]]
					then
						for aa in A R D N C E Q G H I L K M F P S T W Y V
						do
							echo ">"$xmut > ${file%%.fasta}_${rescon}${respos}${aa}.fasta
							sed -r 's/(.{'$position'})(.)/\1'$aa'/' <<< $seqX >> ${file%%.fasta}_${rescon}${respos}${aa}.fasta
						done
					fi
				done
				rm $file
				echo "All sequences from X mutations available."
			fi
		done
	done < ${dataset}-res20_Xmut.log
fi

# Move FASTA to folder
rundate=$(date "+%Y%m%d")
mkdir -p ${rundate}_${dataset}_res20
mv *.fasta ${rundate}_${dataset}_res20

sed '1d' ${dataset}-res15 | while read sequence
do
	seqID=$(awk '{print $1}' <<< $sequence)
	echo ">"$seqID > ${seqID}.fasta
	checkseqID=$(ls | grep $seqID | wc -l)
	
	seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
	# consensus B subtype (https://hivdb.stanford.edu/pages/documentPage/consensus_amino_acid_sequences.html)
	
	for columnno in {10..108}
	do
		aminoacid=$(awk '{print $'$columnno'}' <<< $sequence)
		position=$(expr $columnno - 10)
		
		rescon=${seqcons:$position:1}
		respos=$(expr $position + 1)
		
		checkaa=$(wc -m <<< $aminoacid) # in dataset, if there's a mixture the $(wc -m) is ge 3, else is 2
		if [[ -n $checkseqID && $checkaa -eq 2 ]]
		then
			for file in $(ls | grep $seqID)
			do
				if [[ $aminoacid == '-' ]] ; then printf "%s" $rescon >> $file
				elif [[ $aminoacid == '.' || $aminoacid == '~' ]] ; then
					echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}d.fasta
					printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}d.fasta
				else printf "%s" $aminoacid >> $file ; fi
			done
		elif [[ -n $checkseqID && $checkaa -ge 3 ]]
		then
			for file in $(ls | grep $seqID)
			do
				aa=$(expr $checkaa - 2) # discount for (newline + count from zero)
				for resmix in $(seq 0 $aa)
				do
					resaa=${aminoacid:$resmix:1} # keep mixtures with consensus aminoacid, because of further mutations
					if [[ -n $(grep '*' <<< $resaa) ]] ; then
                                                echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}_STOP.fasta
	                                        printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}_STOP.fasta
					elif [[ -n $(grep '#' <<< $resaa) ]] ; then 
	        				echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
		   				printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
			                        printf "%s" $aminoacid >> ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
                                        else
						echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
						printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
						printf "%s" $resaa >> ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
					fi
				done
				rm $file
			done
		fi
	done
done
echo "All highly resistant sequences greater than 15 are separated in FASTA files."

# Check for X residue
sed '1d' ${dataset}-res15 | while read sequence
do
	seqID=$(awk '{print $1}' <<< $sequence)
	for file in $(ls | grep $seqID)
	do
		checkX=$(sed '1d' $file | grep X)
		if [[ -n $checkX ]] ; then echo $seqID >> ${dataset}-res15_Xmut.tmp ; fi
	done
done
if [[ -f ${dataset}-res15_Xmut.tmp ]]
then
	cat PI_dataset-res15_Xmut.tmp | sort | uniq > PI_dataset-res15_Xmut.log
	rm PI_dataset-res15_Xmut.tmp
	echo "The following sequences have at least one mutation for any possible amino acid:"
	echo
	cat ${dataset}-res15_Xmut.log
	echo
	echo "The IDs are saved in the "${dataset}-res15_Xmut.log" file."

# Separate various X mutations
	while read xmut
	do
		for file in $(ls | grep $xmut)
		do
			seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
			seqX=$(sed '1d' $file | grep X)
			seqNX=$(sed '1d' $file | grep -o X | wc -l)
			if [[ $seqNX -eq 1 ]]
			then
				for position in {0..98}
				do
					resseq=${seqX:$position:1}
					rescon=${seqcons:$position:1}
					respos=$(expr $position + 1)

					## Amino acids
					# "A" = "ALA"
					# "R" = "ARG"
					# "D" = "ASP"
					# "N" = "ASN"
					# "C" = "CYS"
					# "E" = "GLU"
					# "Q" = "GLN"
					# "G" = "GLY"
					# "H" = "HIS"
					# "I" = "ILE"
					# "L" = "LEU"
					# "K" = "LYS"
					# "M" = "MET"
					# "F" = "PHE"
					# "P" = "PRO"
					# "S" = "SER"
					# "T" = "THR"
					# "W" = "TRP"
					# "Y" = "TYR"
					# "V" = "VAL"

					if [[ $resseq == "X" ]]
					then
						for aa in A R D N C E Q G H I L K M F P S T W Y V
						do
							echo ">"$xmut > ${file%%.fasta}_${rescon}${respos}${aa}.fasta
							sed -r 's/(.{'$position'})(.)/\1'$aa'/' <<< $seqX >> ${file%%.fasta}_${rescon}${respos}${aa}.fasta
						done
					fi
				done
				rm $file
				echo "All sequences from X mutations available."
			fi
		done
	done < ${dataset}-res15_Xmut.log
fi

# Move FASTA to folder
rundate=$(date "+%Y%m%d")
mkdir -p ${rundate}_${dataset}_res15
mv *.fasta ${rundate}_${dataset}_res15

sed '1d' ${dataset}-res10 | while read sequence
do
	seqID=$(awk '{print $1}' <<< $sequence)
	echo ">"$seqID > ${seqID}.fasta
	checkseqID=$(ls | grep $seqID | wc -l)
	
	seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
	# consensus B subtype (https://hivdb.stanford.edu/pages/documentPage/consensus_amino_acid_sequences.html)
	
	for columnno in {10..108}
	do
		aminoacid=$(awk '{print $'$columnno'}' <<< $sequence)
		position=$(expr $columnno - 10)
		
		rescon=${seqcons:$position:1}
		respos=$(expr $position + 1)
		
		checkaa=$(wc -m <<< $aminoacid) # in dataset, if there's a mixture the $(wc -m) is ge 3, else is 2
		if [[ -n $checkseqID && $checkaa -eq 2 ]]
		then
			for file in $(ls | grep $seqID)
			do
				if [[ $aminoacid == '-' ]] ; then printf "%s" $rescon >> $file
				elif [[ $aminoacid == '.' || $aminoacid == '~' ]] ; then
					echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}d.fasta
					printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}d.fasta
				else printf "%s" $aminoacid >> $file ; fi
			done
		elif [[ -n $checkseqID && $checkaa -ge 3 ]]
		then
			for file in $(ls | grep $seqID)
			do
				aa=$(expr $checkaa - 2) # discount for (newline + count from zero)
				for resmix in $(seq 0 $aa)
				do
					resaa=${aminoacid:$resmix:1} # keep mixtures with consensus aminoacid, because of further mutations
					if [[ -n $(grep '*' <<< $resaa) ]] ; then
                                                echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}_STOP.fasta
	                                        printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}_STOP.fasta
					elif [[ -n $(grep '#' <<< $resaa) ]] ; then 
	        				echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
		   				printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
			                        printf "%s" $aminoacid >> ${file%%.fasta}_${rescon}${respos}${aminoacid}.fasta
                                        else
						echo ">"$seqID > ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
						printf "%s" $(sed '1d' $file) >> ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
						printf "%s" $resaa >> ${file%%.fasta}_${rescon}${respos}${resaa}.fasta
					fi
				done
				rm $file
			done
		fi
	done
done
echo "All highly resistant sequences greater than 10 are separated in FASTA files."

# Check for X residue
sed '1d' ${dataset}-res10 | while read sequence
do
	seqID=$(awk '{print $1}' <<< $sequence)
	for file in $(ls | grep $seqID)
	do
		checkX=$(sed '1d' $file | grep X)
		if [[ -n $checkX ]] ; then echo $seqID >> ${dataset}-res10_Xmut.tmp ; fi
	done
done
if [[ -f ${dataset}-res10_Xmut.tmp ]]
then
	cat PI_dataset-res10_Xmut.tmp | sort | uniq > PI_dataset-res10_Xmut.log
	rm PI_dataset-res10_Xmut.tmp
	echo "The following sequences have at least one mutation for any possible amino acid:"
	echo
	cat ${dataset}-res10_Xmut.log
	echo
	echo "The IDs are saved in the "${dataset}-res10_Xmut.log" file."

# Separate various X mutations
	while read xmut
	do
		for file in $(ls | grep $xmut)
		do
			seqcons=PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
			seqX=$(sed '1d' $file | grep X)
			seqNX=$(sed '1d' $file | grep -o X | wc -l)
			if [[ $seqNX -eq 1 ]]
			then
				for position in {0..98}
				do
					resseq=${seqX:$position:1}
					rescon=${seqcons:$position:1}
					respos=$(expr $position + 1)

					## Amino acids
					# "A" = "ALA"
					# "R" = "ARG"
					# "D" = "ASP"
					# "N" = "ASN"
					# "C" = "CYS"
					# "E" = "GLU"
					# "Q" = "GLN"
					# "G" = "GLY"
					# "H" = "HIS"
					# "I" = "ILE"
					# "L" = "LEU"
					# "K" = "LYS"
					# "M" = "MET"
					# "F" = "PHE"
					# "P" = "PRO"
					# "S" = "SER"
					# "T" = "THR"
					# "W" = "TRP"
					# "Y" = "TYR"
					# "V" = "VAL"

					if [[ $resseq == "X" ]]
					then
						for aa in A R D N C E Q G H I L K M F P S T W Y V
						do
							echo ">"$xmut > ${file%%.fasta}_${rescon}${respos}${aa}.fasta
							sed -r 's/(.{'$position'})(.)/\1'$aa'/' <<< $seqX >> ${file%%.fasta}_${rescon}${respos}${aa}.fasta
						done
					fi
				done
				rm $file
				echo "All sequences from X mutations available."
			fi
		done
	done < ${dataset}-res10_Xmut.log
fi

# Move FASTA to folder
rundate=$(date "+%Y%m%d")
mkdir -p ${rundate}_${dataset}_res10
mv *.fasta ${rundate}_${dataset}_res10

echo -e "\aIt's finished!"
exit 0
