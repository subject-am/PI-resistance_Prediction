#!/bin/bash

#
#  pattern-HIVp.sh
#
#     Usage:   sh pattern-HIVp.sh -h
#              bash pattern-HIVp.sh 1NH0.pdb pattern01
#
#  Loops the script substitution-HIVp.sh through all the mutations found in
#  the pattern file, to substitute them in the PDB structure found in
#  the present directory and mutates for the indicated substitution in the.
#
#     Example:   sh pattern-HIVp.sh 1NH0.pdb pattern01
#                    keeps intermediary files of the mutation pattern applied to 1NH0.pdb
#                    (assumes 'pattern01' is a list of mutation nomenclatures in separate newlines)
#

# Help Menu
if [[ "$1" == "-h" ]]
then
	echo
	echo "    This script automates the use of 'mutate_model.py' for applying mutation patterns"
	echo "    on an HIV1p structure found in the present working directory and eliminates the"
	echo "    'HETATM' lines from the mutated PDB file. All new files are saved into a new folder."
	echo
	echo "    	pattern-HIVp.sh [HELP]"
	echo "    	pattern-HIVp.sh [ARGS1-3]"
	echo
	echo "    [HELP]	-h	shows this help text"
	echo "    [ARGS1,2] <PDB filename> <mutation pattern exact filename>"
	echo
	echo "    Examples:"
	echo "		./pattern-HIVp.sh 1NH0 pattern01"
	echo "		sh pattern-HIVp.sh 1NH0.pdb pattern01"
	echo
    exit 0
fi

# Get info
pdbf=$1
patf=$2
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )" #source directory of the script

# New directory
patfoo=${pdbf%%.pdb}_${patf%.*}
mkdir -p $patfoo && cd $patfoo

# Create mutation over mutation
while read mutation
do
	sq=$(ls)
	if [[ !  -z  $sq  ]]
	then # it's not empty
		opdb=$npdb
	else # it's empty
		cp ../$pdbf .
		opdb=$pdbf
	fi
	sh ${DIR}/substitution-HIVp.sh -m $mutation $opdb n
	mv mut_${opdb%%.pdb}_${mutation}/${opdb%%.pdb}_${mutation}.pdb .
	rm -r mut_${opdb%%.pdb}_${mutation}/
	npdb="${opdb%%.pdb}_${mutation}.pdb"
done < ../$patf
cd ..
echo -e "\aAll done, "$USER"! A new PDB with a mutation pattern was created in $patfoo.\a"
