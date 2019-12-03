# PI-resistance_Prediction

## Reference Structure
1NH0_ref.pdb was created from 1NH0.pdb (https://files.rcsb.org/download/1NH0.pdb) as follows below.

After downloading the 1nh0.fasta (http://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=1NH0&compressionType=uncompressed), the pattern file for the reference was created from:

        sh pattern2consensus.sh
        
If a reference has no mutations from the consensus subtype B, then:

        echo P1P > ref_pattern

As such, they are converted to the consensus reference template, when subjected to `pattern-HIVp.sh`

which automates the use of `mutate-model.py` (https://salilab.org/modeller/wiki/Mutate%20model).

        sh pattern-HIVp.sh 1NH0.pdb ref_pattern
        mv 1NH0_ref_pattern/$(ls 1NH0_ref_pattern/ --sort=time | head -n 1) 1NH0_ref.pdb
        rm -r 1NH0_ref_pattern/

This removes *altlocs* and HETATM. With the finished model of the reference - and a set *random_seed* - MIF comparison between the models obtained will be possible. Otherwise the grid is not the same.

## Scripts
MutModels.sh

        Will model all the files, performing the necessary mutations to the consensus reference structure.
        
mutate-model.py

        Retrieved from Modeller Wiki (https://salilab.org/modeller/wiki/Mutate%20model).
        Performs single-point substitution mutations on a PDB model.
        
substitution-HIVp.sh

        Automates the use of 'mutate_model.py' for the two chains of an HIV-1 protease
	structure in the present working directory and eliminates the 'HETATM' lines"
	from the mutated PDB file. All new files are saved into a new folder.
        
pattern-HIVp.sh

        Loops the script substitution-HIVp.sh through all the mutations found in
        the pattern file (mutation nomenclature, different lines), to substitute them
        in the PDB structure found in the present directory.
        
pattern2consensus.sh

        Creates the mutation pattern files from consensus B HIV-1 protease.
	
ref_pattern.sh

	Runs the commands specified in the section above to create the reference structure.
