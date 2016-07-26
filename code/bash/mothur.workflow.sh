# Initiate mothur
mothur

# Generate the paired file with the actual sample names
make.file(inputdir=data/WAVES/dataRaw/, outputdir=data/WAVES/mothur, type=fastq)

# Rename the paired file so that it means something for your dataset.
system(mv data/WAVES/mothur/fileList.paired.file data/WAVES/mothur/WAVES.file)

# Make contigs from the sequences generated from the forward and reverse primers.
make.contigs(file=data/WAVES/mothur/WAVES.file, outputdir=data/WAVES/mothur/, processors=2)

# Initial filtering of sequences. Eliminate sequences longer than 275bp, those with homopolymers over 8bp, and those with any ambiguous basepairs.
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=8)

# Retrieve unique sequences
unique.seqs(fasta=current)

# Generate a count table that lists the copy number of each unique sequence)
count.seqs(name=current, group=current)

# Align our sequences to the alignment from the Silva database previously generated.
align.seqs(fasta=current, reference=data/WAVES/mothur/silva.v4.reference.fasta)

# Trim the alignment. Must start before or at 1968, end at or after 11550, and have a max homopolymer length of 8
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

# Classify the sequences using the RDP training database
classify.seqs(fasta=current, count=current, reference=data/references/gg_13_8_99.fasta, taxonomy=data/references/gg_13_8_99.gg.tax, cutoff=80)

# Remove the sequences that are marked at chloroplast, mitochondria, unknown, and eukaryotes.
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Eukaryota)

# Cluster the OTUs, split into their orders first to make it faster (taxlevel=4 signifies Order)
cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.15)

# Make the shared table to see how many of the OTUs are in each sample
make.shared(list=current, count=current, label=0.03)

# Generate a consensus classfication for each OTU
classify.otu(list=current, count=current, taxonomy=current, label=0.03)

# Count em!
count.groups(shared=current)

# Subsample your sequences to the lowest count
# Here, B1014 has the lowest number of sequences, which is 4421
sub.sample(shared=current, size=4421)

# Generate a heatmap of which OTUs are present.
# The subsequent file can be read in Safari?
heatmap.bin(shared=current, scale=log2, numotu=50)

# Generate distance matrix
asd

# quit()
