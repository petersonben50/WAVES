# Directory: /Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES
# Most of the generated files will be stored on BDPLabHD, in
# dataEdited/16STaxAss

# The final output file WAVES.taxonomy is copied to the dataEdited folder in
# the WAVES project directory on my computer

# Remove old mothur log files before starting
rm mothur.*.logfile

# Blast OTUs against FW blast database
blastn -query ./dataEdited/WAVES.final.fasta -task megablast -db /Volumes/BDPLabHD/data/referenceDB/FWonly_11Feb2016_1452_ready.fasta.db -out /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast -outfmt 11 -max_target_seqs 5

# Reformat the blast results
blast_formatter -archive /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast.table

# Use Robin's script to calculate the full length pident of each blast result.
Rscript code/R/calc_full_length_pident.R /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast.table /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast.table.modified

# Use Robin's script to pull out sequence IDs that have greater than 98%
# identity to FW database entry
Rscript code/R/filter_seqIDs_by_pident.R /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast.table.modified /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.above.98 98 TRUE

# Use Robin's script to pull out sequence IDs that have less than 98%
# identity to FW database entry
Rscript code/R/filter_seqIDs_by_pident.R /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast.table.modified /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.below.98 98 FALSE

# Generate a sanity check plot that can be seen in R
# First make a directory.
mkdir /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/plots
# Run Robin's R script to generate the plots.
RScript code/R/plot_blast_hit_stats.R /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast.table.modified 98 /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/plots

# Python code to pull out OTUs that hit nothing on the blast
python code/python/find_seqIDs_blast_removed.py ./dataEdited/WAVES.final.fasta /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/OTU.custom.blast.table.modified /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.missing

# Combine the list of sequence IDs that didn't hit anything on the blast with
# the sequence IDs that had a pident below 98.
cat /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.below.98 /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.missing > /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.below.98.all

# Create a fasta file with the sequences that correspond to each group (above
# or below 98%). This also uses a python code written by Robin.
python code/python/create_fastas_given_seqIDs.py /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.above.98 dataEdited/WAVES.final.fasta /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/otus.above.98.fasta
python code/python/create_fastas_given_seqIDs.py /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.below.98.all dataEdited/WAVES.final.fasta /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/otus.below.98.fasta

# Classify each group of OTUs using mothur with the two different databases.
# Sequences with pident above 98 are classified using the FW database.
# Sequences with pident below 98 are classified using the Gg database.
mothur "#classify.seqs(fasta=/Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/otus.above.98.fasta, template=/Volumes/BDPLabHD/data/referenceDB/FW_Robin.fasta, taxonomy=/Volumes/BDPLabHD/data/referenceDB/FW_Robin.taxonomy, method=wang, probs=T, processors=2)"
mothur "#classify.seqs(fasta=/Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/otus.below.98.fasta, template=/Volumes/BDPLabHD/data/referenceDB/Gg_Robin.fasta, taxonomy=/Volumes/BDPLabHD/data/referenceDB/Gg_Robin.taxonomy, method=wang, probs=T, processors=2)"

# Concatonate the two separate classifications together.
cat /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/otus.above.98.FW_Robin.wang.taxonomy /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/otus.below.98.Gg_Robin.wang.taxonomy > /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy

# Reformat the OTU classification files so that it is delimited by semicolons,
# rather than spaces or tabs.
sed 's/[[:blank:]]/\;/' </Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy >/Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy.reformatted
mv /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy.reformatted /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy

#sed 's/[-]/\_/' </Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy >/Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy.reformatted
#mv /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy.reformatted /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy

# Move a copy of this OTU classification file to dataEdited folder in the WAVES
# directory on my computer.
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/WAVES.taxonomy ./dataEdited/WAVES.taxonomy
