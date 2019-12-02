# PI-resistance_Prediction
## Datasets
PI_dataset

        The complete dataset retrieved from HIVdb from Stanford (last accessed on 20190220 ; file: https://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt) sorted and edited to be read by script.
        Fold-resistance to PI, according to in vitro assay.


Separate_sets.sh

        Reads the PI_dataset and separates the wanted sets, turning them to FASTA files.


PI_dataset-res10

        Filtered sequences that have more than 10.0 fold-resistance to all PI are considered
        to be highly resistance in accordance with https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1622926/


PI_dataset-res15

        Filtered sequences that have more than 15.0 fold-resistance to all PI are considered
        to be highly resistance in accordance with https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1622926/


PI_dataset-res20

        Filtered sequences that have more than 20.0 fold-resistance to all PI are considered
        to be highly resistance in accordance with https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1622926/


PI_dataset-susceptible

        Filtered sequences that have less than 3.0 fold-resistance to all PI are considered
        to be susceptible in accordance with https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1622926/


Resistant_noOut

        Filtered sequences that have more than 10.0 fold-resistance to all PI are considered
        to be highly resistance in accordance with https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1622926/
        without the extreme outliers found through Tukey's method.


Susceptible_noOut

        Filtered sequences that have less than 3.0 fold-resistance to all PI are considered
        to be highly resistance in accordance with https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1622926/
        without the extreme outliers found through Tukey's method.
