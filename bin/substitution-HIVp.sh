#!/bin/bash

#how to use
#       sh substitution-HIVp.sh <PDB filename> <residue position> <residue mutation> <keep intermediary files? yes or no?>
#       sh substitution-HIVp.sh 3DJK 84 VAL n

#help menu
if [[ "$1" == "-h" ]]
then
	echo
	echo "    This script automates the use of 'mutate_model.py' for the two chains of an HIV1p"
	echo "    structure in the present working directory and eliminates the 'HETATM' lines"
	echo "    from the mutated PDB file. All new files are saved into a new folder."
	echo
	echo "    	substitution-HIVp.sh [HELP]"
	echo "    	substitution-HIVp.sh [MUTN] [ARG1]...[ARG4]"
	echo "    	substitution-HIVp.sh [ARGS1-4]"
	echo
	echo "    [HELP]	-h	shows this help text"
	echo "    [MUTN]	-m	uses mutation nomenclature"
	echo "    [ARGS1-4] <PDB filename> <residue position> <mutated residue> <keep intermediary files? yes or no?>"
	echo
	echo "    Examples:"
	echo "		./substitution-HIVp.sh 3DJK 84 VAL y"
	echo "		sh substitution-HIVp.sh 3DJK.pdb 84 val YES"
	echo "		sh substitution-HIVp.sh -m I84V 3DJK.pdb Yes"
	echo
    exit 0
fi

#get info
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" #source directory of the script

if [[ "$1" == "-m" ]]
then
	mn=$3
	mut=${2^^}
	rc=${mut: -1}
	res=${mut%${rc}}
	rp=${res:1}
	of=$4 # works with: Y y YES Yes yes N n NO No no
	if [[ $rc == "A" ]]; then rn="ALA"
	elif [[ $rc == "R" ]]; then rn="ARG"
	elif [[ $rc == "D" ]]; then rn="ASP"
	elif [[ $rc == "N" ]]; then rn="ASN"
	elif [[ $rc == "C" ]]; then rn="CYS"
	elif [[ $rc == "E" ]]; then rn="GLU"
	elif [[ $rc == "Q" ]]; then rn="GLN"
	elif [[ $rc == "G" ]]; then rn="GLY"
	elif [[ $rc == "H" ]]; then rn="HIS"
	elif [[ $rc == "I" ]]; then rn="ILE"
	elif [[ $rc == "L" ]]; then rn="LEU"
	elif [[ $rc == "K" ]]; then rn="LYS"
	elif [[ $rc == "M" ]]; then rn="MET"
	elif [[ $rc == "F" ]]; then rn="PHE"
	elif [[ $rc == "P" ]]; then rn="PRO"
	elif [[ $rc == "S" ]]; then rn="SER"
	elif [[ $rc == "T" ]]; then rn="THR"
	elif [[ $rc == "W" ]]; then rn="TRP"
	elif [[ $rc == "Y" ]]; then rn="TYR"
	elif [[ $rc == "V" ]]; then rn="VAL"
	fi
else
	mn=$1
	rp=$2
	rn=$3
	of=$4 # works with: Y y YES Yes yes N n NO No no
	if [[ ${rn^^} == "ALA" ]]; then rc="A"
	elif [[ ${rn^^} == "ARG" ]]; then rc="R"
	elif [[ ${rn^^} == "ASP" ]]; then rc="D"
	elif [[ ${rn^^} == "ASN" ]]; then rc="N"
	elif [[ ${rn^^} == "CYS" ]]; then rc="C"
	elif [[ ${rn^^} == "GLU" ]]; then rc="E"
	elif [[ ${rn^^} == "GLN" ]]; then rc="Q"
	elif [[ ${rn^^} == "GLY" ]]; then rc="G"
	elif [[ ${rn^^} == "HIS" ]]; then rc="H"
	elif [[ ${rn^^} == "ILE" ]]; then rc="I"
	elif [[ ${rn^^} == "LEU" ]]; then rc="L"
	elif [[ ${rn^^} == "LYS" ]]; then rc="K"
	elif [[ ${rn^^} == "MET" ]]; then rc="M"
	elif [[ ${rn^^} == "PHE" ]]; then rc="F"
	elif [[ ${rn^^} == "PRO" ]]; then rc="P"
	elif [[ ${rn^^} == "SER" ]]; then rc="S"
	elif [[ ${rn^^} == "THR" ]]; then rc="T"
	elif [[ ${rn^^} == "TRP" ]]; then rc="W"
	elif [[ ${rn^^} == "TYR" ]]; then rc="Y"
	elif [[ ${rn^^} == "VAL" ]]; then rc="V"
	fi
fi

oldresidue=$(grep "^ATOM" ${mn%.pdb}.pdb | awk '$6 == '$rp'' | awk '{print $4}' | sed -n '1{p;q}')
oldres=${oldresidue: -3}
if [[ $oldres == "ALA" ]]; then oldr="A"
elif [[ ${oldres} == "ARG" ]]; then oldr="R"
elif [[ ${oldres} == "ASP" ]]; then oldr="D"
elif [[ ${oldres} == "ASN" ]]; then oldr="N"
elif [[ ${oldres} == "CYS" ]]; then oldr="C"
elif [[ ${oldres} == "GLU" ]]; then oldr="E"
elif [[ ${oldres} == "GLN" ]]; then oldr="Q"
elif [[ ${oldres} == "GLY" ]]; then oldr="G"
elif [[ ${oldres} == "HIS" ]]; then oldr="H"
elif [[ ${oldres} == "ILE" ]]; then oldr="I"
elif [[ ${oldres} == "LEU" ]]; then oldr="L"
elif [[ ${oldres} == "LYS" ]]; then oldr="K"
elif [[ ${oldres} == "MET" ]]; then oldr="M"
elif [[ ${oldres} == "PHE" ]]; then oldr="F"
elif [[ ${oldres} == "PRO" ]]; then oldr="P"
elif [[ ${oldres} == "SER" ]]; then oldr="S"
elif [[ ${oldres} == "THR" ]]; then oldr="T"
elif [[ ${oldres} == "TRP" ]]; then oldr="W"
elif [[ ${oldres} == "TYR" ]]; then oldr="Y"
elif [[ ${oldres} == "VAL" ]]; then oldr="V"
fi

#change the python path if necessary
#python='/home/'$USER'/anaconda3/bin/python3.6'
python=$(which python)

####
mkdir -p mut_${mn%.pdb}_${oldr}${rp}$rc
echo "Mutating chain A..."
$python ${DIR}/mutate_model.py ${mn%.pdb} $rp ${rn^^} A > ${mn%.pdb}_A.log
mna=${mn%.pdb}${rn^^}$rp
grep -w "^ATOM" ${mn%.pdb}.pdb > atom.tmp
faa=$(grep -m 1 " B " atom.tmp | awk '{print $6}')
if [ "$faa" == "1" ]
then
	nrp=$rp
elif [ "$faa" == "100" ]
then
	nrp=$(expr 99 + $rp)
else
	nrp=$(expr 100 + $rp)
fi
rm *.tmp
echo "... and chain B..."
$python ${DIR}/mutate_model.py $mna $nrp ${rn^^} B > ${mn%.pdb}_B.log
mnb=$mna${rn^^}$nrp
echo "Deleting ligand..."
sed '/^HETATM/ d' ${mnb}.pdb > ${mn%.pdb}_${oldr}${rp}$rc.pdb
echo "Doing some cleaning..."
if [[ ! $of =~ ^[Yy]$ ]]
then
    rm ${mna}.pdb ${mnb}.pdb
        mv ${mn%.pdb}_* mut_${mn%.pdb}_${oldr}${rp}$rc/
else
        mv ${mna}.pdb ${mnb}.pdb ${mn%.pdb}_* mut_${mn%.pdb}_${oldr}${rp}$rc/
fi
echo -e "\a\a\aIt's done, "$USER"!\nCheck your new "${mn%.pdb}"_"${oldr}${rp}$rc".pdb file.\n"
