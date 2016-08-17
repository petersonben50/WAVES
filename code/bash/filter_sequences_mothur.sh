# This script is to line up the paired end sequences and
# filter out the low quality sequences. It also removes
# chimeras. This is all run on the files in BDPLabHD.

# This script requires a reference database generated from the Silva database.
# The database is a fasta file that contains an alignment of the 16S rRNA segment
# from all the genomes in the Silva database. This file is created in the _____
# script.

# Directory: /Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES
# To run: mothur code/bash/filter_sequences_mothur.sh

# Generate the paired file with the actual sample names
make.file(inputdir=/Volumes/BDPLabHD/data/WAVES/dataEdited/, outputdir=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur, type=gz)

# Rename the paired file so the prefix is always WAVES.
system(mv /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/fileList.paired.file /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.file)

# Make contigs from the sequences generated from the forward and reverse primers.
make.contigs(file=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.file, outputdir=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/, processors=2)

# Initial filtering of sequences. Eliminate sequences longer than 275bp, those
# with homopolymers over 8bp, and those with any ambiguous basepairs.
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=8)

# Retrieve unique sequences
unique.seqs(fasta=current)

# Generate a count table that lists the copy number of each unique sequence
count.seqs(name=current, group=current)

# Align our sequences to the alignment from the Silva database previously generated.
align.seqs(fasta=current, reference=/Volumes/BDPLabHD/data/referenceDB/silva.v4.reference.fasta)

# Trim the alignment. Must start before or at 1968, end at or after 11550,
# and have a max homopolymer length of 8
screen.seqs(fasta=current, count=current, summary=current, start=1968, end=11550, maxhomop=8)

# Filter out the sequences to cut out ones with a period in them
filter.seqs(fasta=current, vertical=T, trump=.)

# Remake the unique sequence list with the filtered sequences
unique.seqs(fasta=current, count=current)

# Precluster the data to denoise it. Speeds up computation
pre.cluster(fasta=current, count=current, diffs=2)

# Remove chimeras
chimera.uchime(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)

# Subsample the shared file so that each sample has the same number of sequences
# in it. This essentially evens out the sampling effort.
# Here, B1014 has the lowest number of sequences, which is 4426
list.seqs(name=current)
get.seqs(group=current, accnos=current)
#sub.sample(fasta=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=current, persample=T,  size=4426)

# Cluster the OTUs.
# First, we need to calculate the distance between sequences
dist.seqs(fasta=current, cutoff=0.20)
# The sequences are then clustered based on that distance matrix
cluster(column=current, count=current)

# Get a representative sequence for each OTU that can then be used to classify
# the OTU. Otherwise, this gets really messy.
get.oturep(column=current, name=current, list=current, fasta=current, label=0.03)

# Make the shared table to see how many of the OTUs are in each sample
#make.shared(list=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.list, count=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, label=0.03)
make.shared(list=current, count=current, label=0.03)
# Get current
#Current files saved by mothur:
#accnos=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos
#column=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.dist
#fasta=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
#group=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.contigs.good.groups
#list=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.list
#name=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.unique.rep.names
#qfile=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.qual
#shared=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared
#count=co
#processors=2
#file=/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/fileList.paired.file
#
#Current input directory saved by mothur: /Volumes/BDPLabHD/data/WAVES/dataEdited/
#
#Current output directory saved by mothur: /Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/
#
#Current default directory saved by mothur: /usr/local/bin/
#
#Current working directory: /Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES/
#
#Output File Names:
#/Volumes/BDPLabHD/data/WAVES/dataEdited/mothur/current_files.summary



# Save a file that has the last current files
get.current()
