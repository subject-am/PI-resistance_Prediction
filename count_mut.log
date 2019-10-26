#!/bin/bash

#### FULL ####

echo "FASTA" > fasta.tmp
echo "#MAJOR" > major.tmp
for f in $(ls *.fasta)
do
        echo ${f%%.fasta} >> fasta.tmp
        checkmajor=$(ls ${f%%.fasta}_majorMut)
        if [[ -z $checkmajor ]] ; then Nmajor=0 ; else Nmajor=$(cat ${f%%.fasta}_majorMut | wc -l) ; fi
        echo $Nmajor >> major.tmp
done
paste fasta.tmp major.tmp > majorMut_histogram.log

echo "FASTA" > fasta.tmp
echo "#MINOR" > minor.tmp
for f in $(ls *.fasta)
do
        echo ${f%%.fasta} >> fasta.tmp
        checkminor=$(ls ${f%%.fasta}_minorMut)
        if [[ -z $checkminor ]] ; then Nminor=0 ; else Nminor=$(cat ${f%%.fasta}_minorMut | wc -l) ; fi
        echo $Nminor >> minor.tmp
done
paste fasta.tmp minor.tmp > minorMut_histogram.log

#### BINDING SITE ####
# (https://hivdb.stanford.edu/pages/3DStructures/pr.html)
# 11-09-2019
# Binding site || Cleft residues: 8,23,25-27,29,30,32,47,48,50,82,84

echo "FASTA" > fasta.tmp
echo "#MAJOR" > majorBS.tmp
for f in $(ls *.fasta)
do
        echo ${f%%.fasta} >> fasta.tmp
        checkmajor=$(ls ${f%%.fasta}_majorMutBS)
        if [[ -z $checkmajor ]] ; then Nmajor=0 ; else Nmajor=$(cat ${f%%.fasta}_majorMutBS | wc -l) ; fi
        echo $Nmajor >> majorBS.tmp
done
paste fasta.tmp majorBS.tmp > majorMut_BS_histogram.log

echo "FASTA" > fasta.tmp
echo "#MINOR" > minorBS.tmp
for f in $(ls *.fasta)
do
        echo ${f%%.fasta} >> fasta.tmp
        checkminor=$(ls ${f%%.fasta}_minorMutBS)
        if [[ -z $checkminor ]] ; then Nminor=0 ; else Nminor=$(cat ${f%%.fasta}_minorMut | wc -l) ; fi
        echo $Nminor >> minorBS.tmp
done
paste fasta.tmp minorBS.tmp > minorMut_BS_histogram.log

rm *.tmp
