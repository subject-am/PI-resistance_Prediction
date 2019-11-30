#!/bin/bash

input=$1 # sh separateFASTA.sh <inputfile>

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
	check=$(ls ${line#>}*.fasta | wc -l)
	fileno=$(expr $check + 1)
	outfile=${line#>}_${fileno}.fasta
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < $input
