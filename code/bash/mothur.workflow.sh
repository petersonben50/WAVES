
# Generate the paired file with the actual sample names
make.file(inputdir=data/WAVES.data/data.raw/, outputdir=data/WAVES.data/mothur, type=fastq)

# Rename that file something meaningful
system(mv data/WAVES.data/mothur/fileList.paired.file data/WAVES.data/mothur/WAVES.file)

# Make contigs from the sequences generated from the forward and reverse primers.
make.contigs(file=data/WAVES.data/mothur/WAVES.file, outputdir=data/WAVES.data/mothur/, processors=2)

# Initial filtering of sequences. Eliminate sequences longer than 275bp, those with homopolymers over 8bp, and those with any ambiguous basepairs.
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=8)

# Retrieve unique sequences
unique.seqs(fasta=current)

# Generate a count table that lists the copy number of each unique sequence)
count.seqs(name=current, group=current)

# Align our sequences to the alignment from the Silva database previously generated.
align.seqs(fasta=current, reference=data/WAVES.data/mothur/silva.v4.reference.fasta)

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
classify.seqs(fasta=current, count=current, reference=data/references/trainset14_032015.pds.fasta, taxonomy=data/references/trainset14_032015.pds.tax, cutoff=80)

# Remove the sequences that obviously don't belong.
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Eukaryota)

# Cluster the OTUs, split into their orders first to make it faster (taxlevel=4 signifies Order)
cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.15)

# Make the shared table to see how many of the OTUs are in each sample
make.shared(list=current, count=current, label=0.03)

# Generate a consensus classfication for each OTU
classify.otu(list=current, count=current, taxonomy=current, label=0.03)

# Count em!
count.groups(shared=current)
