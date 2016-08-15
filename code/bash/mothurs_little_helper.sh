# Mothur's Little Helper

# This script is to rename, move, and reformat files from our mothur folder to
# our dataEdited folder for classification by 16STaxAss. The only files we
# need are the shared files and the fasta files with the unique counts.
# We'll also need to convert the shared file to an OTU file.

# First we'll rename and copy the files into the dataEdited file on BDPLabHD
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.subsample.shared /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.shared
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.fasta

# Let's also move these files onto my computer in the dataEdited directory
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.shared ./dataEdited/WAVES.final.shared
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.fasta ./dataEdited/WAVES.final.fasta

# Now reformat the shared file into an OTU table using an R script that I wrote.
RScript code/R/convert_shared_file.R

# Also add this file to the dataEdited file on BDPLabHD
cp dataEdited/WAVES.final.abund /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.abund
