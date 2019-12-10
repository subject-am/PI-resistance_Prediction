#!/bin/bash
#### Usage || Help menu
echo > usage.help
echo "    This script runs the workflow. Use FASTA for the monomer of" >> usage.help
echo "    HIV-1 protease (subtype B) in one line." >> usage.help
echo "    The script needs Python/Modeller, Reduce, Lovoalign, and IsoMIF to work." >> usage.help
echo >> usage.help
echo "    Example of FASTA file:" >> usage.help
echo "    > Consensus HIV-1 PR (subtype B)" >> usage.help
echo "    PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF" >> usage.help
echo >> usage.help
echo "          sh HIV1predict.sh [HELP]" >> usage.help
echo "          sh HIV1predict.sh [IN] <FASTA file> [PDB] [ARG] [OUT] <output file>" >> usage.help
echo >> usage.help
echo "    [HELP]     -h    shows this help text" >> usage.help
echo "    [IN]       -i    input filename (FASTA format; one sequence only)" >> usage.help
echo "    [PDB]      -p    input own PDB file instead of 1NH0_ref; optional" >> usage.help
echo "                     (assumes that the HIV-1 protease monomer is in chain A)" >> usage.help 
echo "    [ARG]      		  ex. for -p: '3DJK' >> usage.help
echo "    [OUT]      -o    output filename (outputs PDB format)" >> usage.help
echo >> usage.help
echo "    Examples:" >> usage.help
echo "          ./HIV1predict.sh -p 3DJK -o output.pdb -i input.fasta" >> usage.help
echo "          sh HIV1predict.sh -i test -p 1NH0_ref -o test" >> usage.help
echo "          sh HIV1predict.sh -i test.fasta -o test_output.pdb" >> usage.help

usage=$(cat usage.help)

#### Get info from options selected
while getopts i:o:ph option
do
  	case "${option,,}" in
                i) inp=${OPTARG%.fasta};;
                o) out=${OPTARG%.pdb};;
		p) pdb=${OPTARG%.pdb};;
                h) echo "$usage" && rm usage.help && exit 0;;
        esac
done
shift $(( OPTIND -1 ))

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" # source directory of the script

spinner=1 ; spinnerloop="/-\|" ; echo -n ' '

#### Folders where you can find the apps' scripts
reducedir=/applic/reduce/reduce-3.23
isomifdir=/applic/IsoMif/IsoMif_150311
lovodir=/applic/lovoalign/lovoalign-16.342/bin

#### Check apps
if [[ -z $reducedir ]] ; then echo "Check Reduce directory" && rm usage.help && exit 1
elif [[ -z $gmxdir ]] ; then echo "Check GROMACS directory" && rm usage.help && exit 1
elif [[ -z $isomifdir ]] ; then echo "Check IsoMIF directory" && rm usage.help && exit 1
else echo "Directories OK."
echo
fi

#### Let's start
if [[ -n $(ls ${inp}.fasta) ]] ; then input=${inp}.fasta ;
else input=${inp} ; fi
output=${out}.pdb
if [[ -z $pdb ]] ; then pdbfile="${DIR}"/1NH0_ref.pdb ; ref=1NH0_ref.pdb
else
	if [[ -n $(ls ${pdb}.pdb) ]] ; then pdbfile=${pdb}.pdb
	elif [[ -n $(ls ${pdb}) ]] ; then pdbfile=${pdb}
	else wget "https://files.rcsb.org/download/${pdb}.pdb" -O ${pdb}.pdb ; pdbfile=${pdb}.pdb ; fi
fi

#### New working directory
currdir=$(pwd)
inpf=${input##*/}
outf=${output##*/}
pdbf=${pdbfile##*/}
workdir=${inpf%%.fasta}
mkdir -p ./$workdir && cp $input ${workdir}/$inpf && cp $pdbfile ${workdir}/$pdbf
cd ${workdir}/

#### Prepare the FASTA input in one line
## Is there more than one FASTA inside file?
check=$(grep ">" $inpf | wc -l)
if [[ $check -gt 1 ]] ; then
        sh "${DIR}"/separateFASTA.sh $inpf
        ls -I $inpf > seqfasta.tmp
	echo "    The input had more than one FASTA. Please run this script again,"
        echo "    using the new separate FASTA files."
        echo
	grep "    _[[:digit:]]*.fasta" seqfasta.tmp
        echo
	rm seqfasta.tmp
        rm usage.help
        exit 0
fi
##

#### Get the FASTA sequence in one line
tail -c1 $inpf | read -r _ || echo >> $inpf
seqfasta=""
sed '1d' $inpf | while read sequence
do
  	seqfasta=$seqfasta$sequence
	echo $seqfasta > seq.tmp
done
sed -i '2,$d' $inpf && cat seq.tmp >> $inpf
#sed -i '2,$d' $inpf && seq99=$(echo $seqfasta | colrm 100) && echo $seq99 >> $inpf # to keep only 99 aminoacids (beware of insertions!!)

#### Prepare PDB_ref
if [[ $pdbf != 1NH0_ref.pdb ]]
then
	ref=${pdbf%%.pdb}_ref.pdb
	mkdir refTMP ; cd refTMP/
	cp ../$pdbf .
	echo ">"${pdbf%.pdb} > ref.fasta
	# Get PDB to FASTA
	awk '/^ATOM/ && $3 == "CA" && $5 == "A" {print $4}' $pdbf > ref.tmp
	while read aa
	do
		if [[ ${#aa} -gt 3 ]] ; then if [[ ${aa:0:1} == "A" ]] ; then echo ${aa:1} >> ref-fasta.tmp ; fi
		else echo $aa >> ref-fasta.tmp ; fi
	done < ref.tmp
	cat ref-fasta.tmp | tr '\n' ' ' | sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g' | sed 's/ //g' >> ref.fasta
	echo >> ref.fasta
	sh "${DIR}"/mutations-consensus.sh >> build_ref.log
	sh "${DIR}"/pattern-HIVp.sh $pdbf ref_pattern >> build_ref.log
	mv *_ref_pattern/$(ls *_ref_pattern/ --sort=time | head -n 1) ../$ref
	cd .. ; rm -r refTMP/
	echo "    Reference model finished."
	echo
fi

#### Build model
echo -n "    Adding mutations...    "
while [ ! -f ${outf%%.pdb}.pdb ] ; do printf "\b${spinnerloop:spinner++%${#spinnerloop}:1}" ; sleep 0.3 ; done &
#while [ ! -f ${outf%%.pdb}.pdb ]; do echo -e -n '|\r\c' ; sleep 0.3 ; echo -e -n '/\r\c' ; sleep 0.3 ; echo -e -n '-\r\c' ; sleep 0.3 ; echo -e -n '\\\r\c' ; sleep 0.3 ; done &
sh "${DIR}"/mutations-consensus.sh >> build_model.log
sh "${DIR}"/pattern-HIVp.sh $ref *_pattern >> build_model.log
mv *ref_*_pattern/$(ls *ref_*_pattern/ --sort=time | head -n 1) ./${outf%%.pdb}.pdb # most recent
rm -r *ref_*_pattern/
echo
echo "    Structural model finished."
echo

#### Prepare reference for IsoMIF
echo -n "    Starting IsoMIF...     "
while [ ! -f isomif_models/*.isomif ] ; do printf "\b${spinnerloop:spinner++%${#spinnerloop}:1}" ; sleep 0.3 ; done &
mkdir isomif_ref && cd isomif_ref
${reducedir}/reduce ../$ref > ${ref%%.pdb}_H.pdb 2>> mif_ref.log	# add H
${isomifdir}/getcleft_linux_x86_64 -p ${ref%%.pdb}_final.pdb -s -t 1 -k 1 -o ${ref%%.pdb}_final &>> mif_ref.log 	# getCleft
${isomifdir}/mif_linux_x86_64 -p ${ref%%.pdb}_final.pdb -g ${ref%%.pdb}_final_sph_1.pdb -t ${ref%%.pdb}_final &>> mif_ref.log  # obtain MIFs for the chosen cleft
cd ..

#### IsoMIF
mkdir isomif_models && cd isomif_models
referencedir=../isomif_ref
${reducedir}/reduce ../${outf%%.pdb}.pdb > ${outf%%.pdb}_H.pdb 2>> mif_models.log # add H
${lovodir}/lovoalign -p1 ${outf%%.pdb}_H.pdb -p2 ${referencedir}/1NH0_ref_final.pdb -o ${outf%%.pdb}_final.pdb &>> mif_models.log	# align model and reference
${isomifdir}/mif_linux_x86_64 -p ${outf%%.pdb}_final.pdb -g ${referencedir}/*_final_sph_1.pdb -t ${outf%%.pdb}_mif &>> mif_models.log
$isomifdir/isomif_linux_x86_64 -p1 ${referencedir}/*_final.mif -p2 ${outf%%.pdb}_mif.mif -s 1 -c 1 -d 1.0 -w -o grid1_d10_ &>> isomif.log
cd ..
echo
echo "    IsoMIF calculated."
echo

#### Resistant? Susceptible?
isomif=$(ls isomif_models/*.isomif)
alltanimoto=$(grep TANI $isomif | tr " " "\n" | grep TANI -n | awk -F ":" '{print $1}')
for coln in $alltanimoto
do
	tcol=$(grep TANI $isomif | awk '{print $'$coln'}')
	if [[ $tcol == "TANI" ]] ; then tanicol=$(( $coln + 1 )) ; tani=$(grep TANI $isomif | awk '{print $'$tanicol'}') ; fi
done
disscoef=$(awk -v n1=1 -v n2=$tani -v OFMT="%.4f" 'BEGIN{print n1-n2}')

if [[ $disscoef > 0.0600 ]] ; then result=resistant ; else result=susceptible ; fi
echo "    Your "$inpf" has a dissimilarity of "$disscoef" compared to"
echo "    the susceptible reference."
echo "    This sequence is likely "$result" to most PIs."
echo

#### Time and Clean-up
times &>> time.tmp
tail -n1 time.tmp | awk '{print $1}'
echo

rm *.tmp
cd $currdir
rm usage.help
exit 0
