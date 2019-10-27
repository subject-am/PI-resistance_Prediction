#!/bin/bash

	# 1NH0_ref.pdb was created from by eliminating the I and J chains (polypeptides that served as inhibitor/ligands)
	# Pattern files for the reference was created from:
sh pattern2consensus.sh
	# if a reference has no mutations from the consensus subtype B, then it will mutate P1P
	# As such, they were reversed to create the consensus reference template.
sh pattern-HIVp.sh 1NH0.pdb ref_pattern
mv 1NH0_ref_pattern/$(ls 1NH0_ref_pattern/ --sort=time | head -n 1) 1NH0_ref.pdb
rm -r 1NH0_ref_pattern/
