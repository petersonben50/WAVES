# Overall workflow

# Move the .gz files.
bash code/bash/mise_en_place.sh

# Code to generate needed databases?

# Align the sequences and filter them using mothur
bash filter_sequences_mothur.sh

# Move the shared and fasta files. Convert the fasta file
# to an OTU table using an R script
bash mothurs_little_helper.sh

# Run the 16STaxAss workflow.
bash code/bash/16STaxAss_workflow.sh
