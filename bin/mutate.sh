#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" #source directory of the script

#### Usage || Help menu
echo > usage.help
echo "    This script helps mutate FASTA sequences. Use aminoacid codes." >> usage.help
echo >> usage.help
echo "          sh mutate.sh [HELP]" >> usage.help
echo "          sh mutate.sh [IN] <FASTA file> [MUTN] [ARG1] .. [MUTN] [ARGN] [OUT] <output file>" >> usage.help
echo >> usage.help
echo "    [HELP]     -h    shows this help text" >> usage.help
echo "    [IN]       -i    input filename (FASTA format; one sequence only)" >> usage.help
echo "    [MUTN]     -s    substitution mutation (repeat for different mutations)" >> usage.help
echo "               -a    insertion mutation (works with only one mutation)" >> usage.help
echo "    [ARGS1..]        ex. for -s: 'D30N' ; ex. for -a: 'E35E_T'" >> usage.help
echo "    [OUT]      -o    output filename (outputs FASTA format)" >> usage.help
echo >> usage.help
echo "    Examples:" >> usage.help
echo "          ./mutate.sh -s D30N -o 2R8N_mut.fasta -i 2R8N.fasta" >> usage.help
echo "          sh mutate.sh -i 2R8N.fasta -s D30N -s I47V -o 2R8N_mut" >> usage.help
echo "          sh mutate.sh -i 2R8N.fasta -s D30N -s I47V -a E35E_T -o 2R8N_mut" >> usage.help

usage=$(cat usage.help)

#### Get info from options selected
while getopts i:s:a:o:h option
do
        case "${option,,}" in
                i) inpf=${OPTARG};;
                s) subs+=("$OPTARG");;
                a) adds=${OPTARG};;
                o) outpf=${OPTARG%.fasta};;
		h) echo "$usage" && rm usage.help && exit 0;;
        esac
done
shift $(( OPTIND -1 ))

#### Prepare the FASTA input in one line
## Is there more than one FASTA inside file?
check=$(grep ">" $inpf | wc -l)
if [[ $check -gt 1 ]] ; then
	sh "${DIR}"/separateFASTA.sh $inpf
	ls -I $inpf > seqfasta.tmp
	echo
	echo "The input had more than one FASTA. Please run this script again,"
	echo "using the new separate FASTA files."
	echo
	grep "_[[:digit:]]*.fasta" seqfasta.tmp
	echo
	rm seqfasta.tmp
	rm usage.help
	exit 0
fi
##

#### Get the FASTA sequence in one line
seqfasta=""
sed '1d' $inpf | while read sequence
do
	seqfasta=$seqfasta$sequence
	echo $seqfasta > seq.tmp
done
sed -i '2,$d' $inpf
cat seq.tmp >> $inpf

#sed -i '2,$d' $inpf && seq99=$(echo $seqfasta | colrm 100) && echo $seq99 >> $inpf # to keep only 99 aminoacids
seqfasta=$(cat seq.tmp)

#### Substitution mutations first, to avoid errors
if [[ -n $subs ]]
then
	for mut in "${subs[@]}"
	do
		rmut=${mut:1}
		resnew=${mut: -1}
		respos=${rmut%%$resnew}
		position=$(expr $respos - 1)
		seqbefore=${seqfasta:0:$position}
		seqafter=${seqfasta:$respos} # 1st character in string is at position 0
		seqfasta=$seqbefore$resnew$seqafter
		echo $seqfasta > seq.tmp
	done
fi
seqfasta=$(cat seq.tmp)

#### Insertion afterwards
if [[ -n $adds ]]
then
	resmut=${adds:1}
	resadd=${adds: -1}
	respos=${resmut%%?_$resadd}
	seqbefore=${seqfasta:0:$respos} # insertion happens at respos+1
	seqafter=${seqfasta:$respos}
	seqfasta=$seqbefore$resadd$seqafter
	echo $seqfasta > seq.tmp
fi	
seqfasta=$(cat seq.tmp)

#### Output
if [ -z $outpf ] ; then outpf=${inpf%.*}_output.fasta ; fi
echo ">${outpf%.fasta}|mutate" > ${outpf%.fasta}.fasta
echo $seqfasta >> ${outpf%.fasta}.fasta
rm *.tmp
rm usage.help
echo -e "\aDone!\n"
exit 0
