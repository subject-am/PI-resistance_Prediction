#!/bin/bash

# (https://hivdb.stanford.edu/dr-summary/resistance-notes/PI/)
# 31-10-2018
# XXXr highest level of resistance to PI
# XXXp potential contraindication in the use of PI

# D30N    NFVr
# V32I    ATV DRVp FPVr IDVp LPVp TPV
# L33F    ATV DRV FPVp LPV NFV TPV
# M46I    ATV FPVp IDVp LPV NFVp TPV
# M46L    ATV FPVp IDVp LPV NFVp TPV
# I47V    ATV DRVp FPVp IDV LPV NFV TPVp
# I47A    DRVp FPVr LPVr TPVp
# G48V    ATVp LPV NFVr SQVr
# G48M    ATVp LPV NFVr SQVr
# I50L    ATVr
# I50V    DRVp FPVr LPVp
# I54V    ATV FPV IDVr LPVp NFVr SQVp TPVp
# I54T    ATV FPV IDVr LPVp NFVr SQVp
# I54A    ATV FPV IDVr LPVp NFVr SQVp TPVp
# I54L    ATV DRVp FPVr IDV LPVp NFVr SQV
# I54M    ATV DRVp FPVr IDV LPVp NFVr SQV TPVp
# L76V    DRVp FPVr IDVr LPVr
# V82A    ATV FPV IDVr LPVr NFVr SQV
# V82F    ATV DRV FPVr IDVr LPVr NFVr
# V82S    ATV FPV IDVr LPVr NFVr
# V82T    ATV FPV IDVr LPVr NFVr SQV TPVr
# V82L    TPVr
# I84V    ATVr DRV FPVr IDVr LPVp NFVr SQVr TPVp
# N88S    ATVr IDV NFVr SQV
# N88D    NFVr
# L90M    ATVp FPVp IDVp LPV NFVr SQVr

# Files to filter
ls *_pattern > pattern.tmp

# Filter resistant mutations
while read patfile
do
	while read patmut
	do
		if [[ $patmut == D30N ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == V32I ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == L33F ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == M46[IL] ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == I47[VA] ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == G48[VM] ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == I50[LV] ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == I54[VTALM] ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == L76V ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == V82[AFSTL] ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == I84V ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == N88[SD] ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		elif [[ $patmut == L90M ]]; then echo $patmut >> ${patfile%%_pattern}_majorMut
		fi
	done < $patfile
done < pattern.tmp
rm pattern.tmp
