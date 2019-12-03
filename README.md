# PI-resistance_Prediction

## Reference Structure
1NH0_ref.pdb was created from 1NH0.pdb (https://files.rcsb.org/download/1NH0.pdb) as follows below.

Pattern files for the reference are created from:

        `sh pattern2consensus.sh`
        
if a reference has no mutations from the consensus subtype B, then it will mutate P1P

As such, they are converted to the consensus reference template, when subjected to `pattern-HIVp.sh`

which automates the use of `mutate-model.py` (https://salilab.org/modeller/wiki/Mutate%20model).

        `sh pattern-HIVp.sh 1NH0.pdb ref_pattern`
        
        `mv 1NH0_ref_pattern/$(ls 1NH0_ref_pattern/ --sort=time | head -n 1) 1NH0_ref.pdb`
        
        `rm -r 1NH0_ref_pattern/`
