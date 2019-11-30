#!/bin/bash
categories=$(ls -d ${PWD}/*/) # absolute path
workdir=${PWD} # absolute path
#categories=$(ls -d *) # relative path
​
for category in $categories
do
	folder=${category%%/}
	cd $folder
	echo "SEQID_PATTERN" > pattern.tmp
	echo "TANI" > tani.tmp
	echo "TANIM" > tanim.tmp
	echo "TANIMW" > tanimw.tmp
	echo "TANINORM" > taninorm.tmp
​
	for isomif in $(ls *.isomif)
	do
		mifpattern=$(grep "REMARK tag2:" $isomif -m 1 | awk '{print $3}')
		fpattern=${mifpattern%%_mif} ; idpattern=${fpattern##3DJK_ref_}
		echo $idpattern >> pattern.tmp
		alltanimoto=$(grep TANI $isomif | tr " " "\n" | grep TANI -n | awk -F ":" '{print $1}')
		for coln in $alltanimoto
		do
			tcol=$(grep TANI $isomif | awk '{print $'$coln'}')
			if [[ $tcol == "TANI" ]] ; then tanicol=$(( $coln + 1 )) ; tani=$(grep TANI $isomif | awk '{print $'$tanicol'}') ; echo $tani >> tani.tmp
			elif [[ $tcol == "TANIM" ]] ; then tanimcol=$(( $coln + 1 )) ; tanim=$(grep TANIM $isomif | awk '{print $'$tanimcol'}') ; echo $tanim >> tanim.tmp
			elif [[ $tcol == "TANIMW" ]] ; then tanimwcol=$(( $coln + 1 )) ; tanimw=$(grep TANIMW $isomif | awk '{print $'$tanimwcol'}') ; echo $tanimw >> tanimw.tmp
			elif [[ $tcol == "TANINORM" ]] ; then taninormcol=$(( $coln + 1 )) ; taninorm=$(grep TANINORM $isomif | awk '{print $'$taninormcol'}') ; echo $taninorm >> taninorm.tmp
			fi
		done
	done
	catname=${folder##*/}
	paste pattern.tmp tani.tmp tanim.tmp tanimw.tmp taninorm.tmp > isomif_${catname}.tsv
	rm *.tmp
	cd $workdir
	#cd .. # relative path
done
