#!/bin/bash

# (https://hivdb.stanford.edu/dr-summary/resistance-notes/PI/)
# 03-09-2019
# XXXc common nonpolymorphic accessory resistance mutations to PI
# XXXa accessory resistance mutations associated with reduced PI susceptibility
# XXXr rare nonpolymorphic PI-selected variants (not well-studied)

# L10F    DRVc FPVc IDVc LPVc NFVc
# L10I    ATVa DRVa FPVa IDVa LPVa NFVa SQVa TPVa
# L10V    ATVa DRVa FPVa IDVa LPVa NFVa SQVa TPVa
# L10R    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# L10Y    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# V11I    DRVa FPVa
# V11L    DRVr FPVr
# K20R    ATVc DRVc FPVc IDVc LPVc NFVc SQVc TPVc
# K20I    NFVa # consensus amino acid in subtype G and CRF02_AG viruses
# K20M    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# K20T    ATVa FPVa IDVa LPVa NFVa SQVa
# K20V    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# L23I    NFVr
# L24I    ATVc FPVc IDVc LPVc SQVc # increases susceptibility to TPV
# L24F    ATVa FPVa IDVa LPVa SQVa
# L24M    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# M36I    ATVa DRVa FPVa IDVa LPVa NFVa SQVa TPVa # increases replication fitness of viruses with PI-resistance mutations # consensus amino acid in most of the non-B subtypes
# K43T    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVa
# M46V    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# G48A    ATVr LPVr NFVr SQVr
# G48S    ATVr LPVr NFVr SQVr
# G48T    ATVr LPVr NFVr SQVr
# G48Q    ATVr LPVr NFVr SQVr
# G48L    ATVr LPVr NFVr SQVr
# F53L    ATVc IDVc LPVc NFVa SQVc
# F53Y    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# I54S    ATVa FPVa IDVa LPVa NFVa SQVa TPVa
# I54T    ATVa FPVa IDVa LPVa NFVa SQVa TPVa
# Q58E    ATVc DRVc FPVc IDVc LPVc NFVc SQVc TPVc
# A71V    ATVc DRVc FPVc IDVc LPVc NFVc SQVc TPVc
# A71T    ATVc DRVc FPVc IDVc LPVc NFVc SQVc TPVc
# A71I    ATVa DRVa FPVa IDVa LPVa NFVa SQVa TPVa
# A71L    ATVa DRVa FPVa IDVa LPVa NFVa SQVa TPVa
# G73S    ATVc DRVa FPVa IDVc LPVa NFVc SQVc TPVa
# G73T    ATVc DRVa FPVa IDVc LPVa NFVc SQVc TPVa
# G73C    ATVc DRVa FPVa IDVc LPVa NFVc SQVc TPVa
# G73A    ATVc DRVa FPVa IDVc LPVa NFVc SQVc TPVa
# G73D    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# G73V    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# T74P    ATVa DRVa FPVa IDVa LPVa NFVa SQVa TPVa
# T74S    ATVa DRVa FPVa IDVa LPVa NFVa SQVa TPVa # consensus amino acid in most non-B subtypes
# V82M    ATVr DRVr FPVr IDVa LPVa NFVr SQVr TPVr # in most subtypes, it requires a 2-bp mutation # in subtype G, it usually requires a single bp mutation (IDV!)
# V82C    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# N83D    ATVc IDVc NFVc SQVc TPVc
# N83S    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# I84A    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr # high-level resistance to PI
# I84C    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr # less marked effect on PI susceptibility
# I85V    ATVc DRVc FPVc IDVc LPVc NFVc SQVc TPVc # minimal, if any, effect on PI susceptibility
# N88T    ATVa DRVr FPVr IDVr LPVr NFVa SQVr TPVr
# N88G    ATVa DRVr FPVr IDVr LPVr NFVa SQVr TPVr
# L89V    DRVc FPVc IDVc LPVc NFVc
# L89T    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr # uncertain phenotypic and clinical significance

# Files to filter
ls *_pattern > pattern.tmp

#### Minor Mutations ####
while read patfile
do
	while read patmut
	do
		if [[ $patmut == L10[FIVRY] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == V11[IL] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == K20[RIMTV] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == L23I ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == L24[IFM] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == M36I ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == K43T ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == M46V ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == G48[ASTQL] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == F53[LY] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == I54[ST] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == Q58E ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		elif [[ $patmut == A71[VTIL] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
                elif [[ $patmut == G73[STCADV] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
                elif [[ $patmut == T74[PS] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
                elif [[ $patmut == V82[MC] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
                elif [[ $patmut == N83[DS] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
                elif [[ $patmut == I84[AC] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
                elif [[ $patmut == I85V ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
                elif [[ $patmut == N88[TG] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
                elif [[ $patmut == L89[VT] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMut
		fi
	done < $patfile
done < pattern.tmp

# L23I    NFVr
# G48A    ATVr LPVr NFVr SQVr
# G48S    ATVr LPVr NFVr SQVr
# G48T    ATVr LPVr NFVr SQVr
# G48Q    ATVr LPVr NFVr SQVr
# G48L    ATVr LPVr NFVr SQVr
# V82M    ATVr DRVr FPVr IDVa LPVa NFVr SQVr TPVr # in most subtypes, it requires a 2-bp mutation # in subtype G, it usually requires a single bp mutation (IDV!)
# V82C    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr
# I84A    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr # high-level resistance to PI
# I84C    ATVr DRVr FPVr IDVr LPVr NFVr SQVr TPVr # less marked effect on PI susceptibility

#### Minor Mutations in the Binding Site ####
while read patfile
do
	while read patmut
	do
		if [[ $patmut == L23I ]]; then echo $patmut >> ${patfile%%_pattern}_minorMutBS
		elif [[ $patmut == G48[ASTQL] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMutBS
                elif [[ $patmut == V82[MC] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMutBS
                elif [[ $patmut == I84[AC] ]]; then echo $patmut >> ${patfile%%_pattern}_minorMutBS
		fi
	done < $patfile
done < pattern.tmp
rm pattern.tmp
