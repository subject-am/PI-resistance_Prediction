#!/bin/bash

for f in $(ls *.fasta)
do
  checkMmut=$(ls ${f%.fasta}* | grep _majorMut)
  checkmmut=$(ls ${f%.fasta}* | grep _minorMut)
  if [[ -n $checkMmut ]] && [[ -n $checkmmut ]]
  then echo ${f%.fasta}_pattern >> majorM_minorM.log      # major and minor HIVdb
  elif [[ -n $checkMmut ]]
  then echo ${f%.fasta}_pattern >> majorM.log     # major HIVdb
  elif [[ -n $checkMmut ]] && [[ -z $checkmmut ]]
  then echo ${f%.fasta}_pattern >> majorM-only.log     # only major HIVdb
  elif [[ -z $checkMmut ]] && [[ -n $checkmmut ]]
  then echo ${f%.fasta}_pattern >> minorM-only.log        # only minor HIVdb
  elif [[ -n $checkmmut ]]
  then echo ${f%.fasta}_pattern >> minorM.log     # minor HIVdb
  elif [[ -z $checkMmut ]] && [[ -z $checkmmut ]]
  then echo ${f%.fasta}_pattern >> SFW.log
  fi
done

for f in $(ls *.fasta)
do
	check=$(ls ${f%%.fasta}_pattern)
	if [[ -z $check ]]
	then echo $f >> consensus.log
	fi
done
