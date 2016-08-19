# This is a clean up file to prep the directory.
# Before this is run, we need a dataEdited and a dataRaw directory.
# The dataRaw directory should contain the raw zipped sequencing files.

# Directory: /Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES
# Run script: bash code/bash/mise_en_place.sh

# Copy the zipped files to the data edited file.
cp $(find /Volumes/BDPLabHD/data/WAVES/dataRaw/ -name *.gz) /Volumes/BDPLabHD/data/WAVES/dataEdited

# Remove any ._ files in the dataEdited file
rm /Volumes/BDPLabHD/data/WAVES/dataEdited/._*

# Remove old mothur log files before starting
# Only needed if I'm rerunning the file.
# Can you tell I hate having unneeded files lying around?
rm mothur.*.logfile

# Prepare the database files that are needed.
makeblastdb -dbtype nucl -in /Volumes/BDPLabHD/data/referenceDB/FWonly_11Feb2016_1452_ready.fasta -input_type fasta -parse_seqids -out /Volumes/BDPLabHD/data/referenceDB/FWonly_11Feb2016_1452_ready.fasta.db
