# Mothur's Little Helper

# This script is to rename, move, and reformat files from our mothur folder to
# our dataEdited folder for classification by 16STaxAss. The only files we
# need are the shared files and the fasta files with the unique counts.
# We'll also need to convert the shared file to an OTU file.

# First we'll rename and copy the files into the dataEdited file on BDPLabHD
# Shared file - To know the counts per sample
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.shared
# Fasta file - Obviously.
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.rep.fasta /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.longname.fasta
# List file - List of sequences that go with each OTU
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.list /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.list
# Name file - For use in the taxonomy summary later
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.rep.names /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.names
# Group file - For use in the taxonomy summary later
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.contigs.good.pick.groups /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.groups

# Let's also move these files onto my computer in the dataEdited directory
# Shared file
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.shared ./dataEdited/WAVES.final.shared
# Fasta file
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.longname.fasta ./dataEdited/WAVES.final.longname.fasta
# List file
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.list ./dataEdited/WAVES.final.list
# Name file - For use in the taxonomy summary later
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.names ./dataEdited/WAVES.final.names
# Group file - For use in the taxonomy summary later
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.final.groups ./dataEdited/WAVES.final.groups

# Reformat the fasta file names so the sequence name is >OtuXXXX
sed 's/.*\(Otu[0-9]*\).*/>\1/' /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.longname.fasta > /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.fasta

# Also move this code onto BDPLabHD
cp /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.fasta ./dataEdited/WAVES.final.fasta

# Now reformat the shared file into an OTU table using an R script that I wrote.
RScript code/R/convert_shared_file.R

# Also add this file to the dataEdited file on BDPLabHD
cp dataEdited/WAVES.final.abund /Volumes/BDPLabHD/data/WAVES/dataEdited/WAVES.final.abund
