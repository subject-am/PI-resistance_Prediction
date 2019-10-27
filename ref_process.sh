#!/bin/bash

for pdb in 1NH0_ref;do
	/applic/reduce/reduce-3.23/reduce -p ${pdb}.pdb > ${pdb}_H.pdb
	/applic/gromacs/gromacs-2016.4/bin/gmx editconf -f ${pdb}_H.pdb -o ${pdb}_final.pdb -center 0 0 0
	/applic/IsoMif/IsoMif_150311/getcleft_linux_x86_64 -p ${pdb}_final.pdb -s -t 1 -k 1 -o ${pdb}_final
	/applic/IsoMif/IsoMif_150311/mif_linux_x86_64 -p ${pdb}_final.pdb -g ${pdb}_final_sph_1.pdb -t ${pdb}_final
done
