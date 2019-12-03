# PI-resistance_Prediction
## Requirements
This project was developed to run on Linux Bash shell.
To run `HIV1predict.sh` as is, the following need to be previously installed:
 - Modeller:     https://salilab.org/modeller/download_installation.html
 - Reduce:       http://kinemage.biochem.duke.edu/software/reduce.php
 - Gromacs:      http://manual.gromacs.org/documentation/
 - LovoAlign:    https://www.ime.unicamp.br/~martinez/lovoalign/software.html
 - IsoMIF:       https://github.com/mtthchrtr/IsoMif

To run this in Windows, consider enabling Linux Subsystem for Linux and install Bash shell on Windows.

## Installing HIV1predict.sh
Clone the repository:

    git clone https://github.com/subject-am/PI-resistance_Prediction.git
    chmod +x PI-resistance_Prediction/bin/HIV1predict.sh

or download the ZIP file of the **master branch**:

    wget https://github.com/subject-am/PI-resistance_Prediction/archive/master.zip
    unzip master.zip
    chmod +x PI-resistance_Prediction-master/bin/HIV1predict.sh
    
Update the directories for Reduce, Gromacs, LovoAlign, and Isomif on `bin/HIV1predict.sh`:

    #### Folders where you can find the apps' scripts
    reducedir=/applic/reduce/reduce-3.23
    gmxdir=/applic/gromacs/gromacs-2016.4/bin 
    isomifdir=/applic/IsoMif/IsoMif_150311
    lovodir=/applic/lovoalign/lovoalign-16.342/bin

and, if necessary, change the path for Python/Modeller on `bin/substitution-HIVp.sh`:

    #change the python path if necessary
    #python='/home/'$USER'/anaconda3/bin/python3.6'
    python=$(which python)
    
## Running HIV1predict.sh
For ease of use, create an alias on .bashrc, .bash_profile, or .profile, or add it to $PATH.

    alias HIVpredict='<install directory>/PI-resistance_Prediction/bin/HIV1predict.sh'

Using an aminoacid FASTA file with just the sequence for the monomer of HIV-1 protease:

    HIV1predict -i <input>.fasta -o <output>.pdb

will run the sequence of aminoacids to search for mutations from consensus B HIV-1 protease, model it, and compare the MIFs of both mutated model and reference, to classify it as susceptible or resistant in _stdout_. Below are the examples of the screen output for the FASTA files in `examples/`:

    [user@examples]$ sh bin/HIV1predict.sh -i 99mut.fasta -o 99mut
    Apps' directories OK.

    Structural model finished.

    Starting IsoMIF...

    IsoMIF calculated.

    Your 99mut.fasta has a dissimilarity of 0.9170 compared to
    the susceptible reference.
    This sequence is likely resistant to most PIs.

    5m46.501s

    [user@examples]$ sh bin/HIV1predict.sh -i consensus.fasta -o consensus
    Apps' directories OK.

    Structural model finished.

    Starting IsoMIF...

    IsoMIF calculated.

    Your consensus.fasta has a dissimilarity of 0.0057 compared to
    the susceptible reference.
    This sequence is likely susceptible to most PIs.

    0m31.679s

If the input file is in multi-FASTA format, the script will separate the sequences in multiple FASTA files and exit. Run again for each of the new FASTA files created. If the sequence of aminoacids is not in one line, the script will create a new file with the sequence in one line and run normally.


Optionally, the user may want to choose their own template:

    HIV1predict -i <input>.fasta -p <template>.pdb -o <output>.pdb
    
For help:

    HIV1predict -h
    
## Repository
### bin/
Contains the 1NH0_ref.pdb (reference structure) and scripts necessary to run `HIV1predict.sh`:
  - `separateFASTA.sh`: separates multi-FASTA format files.
  - `mutations-consensus.sh`: creates the mutation pattern file.
  - `pattern-HIVp.sh`: creates the mutation pattern file; uses `substitution-HIVp.sh` and `mutate_model.py`.

All of these can be run alone. Other two scripts in `bin/` are not necessary, but can be helpful:
  - `mutate.sh`: mutates a FASTA sequence (deals with FASTA as mentioned above).
  - `isomif_summary.sh`: summarises all Tanimotos resulting from IsoMIF.


### examples/
Examples of the folders and files created when running `HIV1predict.sh` for **99mut.fasta** and **consensus.fasta**.
