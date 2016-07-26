#### 2016-07-16

Outline of pipeline

1. Move files into a usable working directory
2. Use FastQC to check the quality of the data.
3. Remove any sequences that do not align properly.
4. Move into mothur. Align the paired-end sequences and filter for quality
5. Classify them using Green Genes.

This will be run in an Amazon instance.
    1. t2.micro
    2. us-east-1c
    3. Connect to my Volume on EC2.

##### Quality Checking

1. Get FastQC and get it so it can be used (need to get java as well)

```
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip

sudo apt-get install unzip
unzip fastqc_v0.11.5.zip
chmod 755 FastQC/fastqc
sudo add-apt-repository ppa:webupd8team/java
sudo apt-get update
sudo apt-get install oracle-java7-installer
```

2. Get access to data

```
sudo mkdir /data
sudo mount /dev/xvdf /data
```

3. Run FastQC on all samples

```
~/FastQC/fastqc /data/WAVES_data/data_raw/*.fastq -o /data/WAVES_data/fastqc_output/
```

4. Use a for loop to unzip all FastQC files

```
for i in /data/WAVES_data/fastqc_output/*.zip
do
  echo "Unzipping $i"
  unzip $i -d /data/WAVES_data/fastqc_output/
done

cat /data/WAVES_data/fastqc_output/*fastqc/summary.txt > /data/WAVES_data/fastqc_output/fastqc_output.txt
```

Look at this later. Let's get into mothur!!!!

5. Need to go into mothur-enabled AMI
Use mothur AMI, m3.large
Named: WAVES_mothur_session
Security: keep the default "ssh" rule, added the following rules by clicking on "add rule" to add everything below
    a. HTTP (no change)
    b. HTTPS (no change)
    c. custom TCP, change "Port Range" to 8787, change "Source" to "Anywhere"
ssh -i ~/Downloads/edamameShellLesson.pem ubuntu@ec2-52-91-15-242.compute-1.amazonaws.com

6. Connect to volume
```
sudo mkdir data
sudo mount /dev/xvdf data
```

7. Set up mothur directories (not sure if the passwords are needed)
```
mkdir data/WAVES_data/mothur
mkdir data/WAVES_data/data_edited
chmod -R 0777 data/WAVES_data/mothur/
chmod -R 0777 data/WAVES_data/data_raw/
chmod 755 data/WAVES_data/data_raw/
chmod 755 data/WAVES_data/mothur/
```

8. Start mothur pipeline
```
mothur
```
Make the make.file

```
make.file(inputdir=data/WAVES_data/data_raw, outputdir=data/WAVES_data/mothur, type=fastq)
```

Rename the file containing information about which one is paired.

```
system(mv data/WAVES_data/mothur/fileList.paired.file data/WAVES_data/mothur/WAVES.file)
```

Make the contigs from the forward and reverse primers.

```
make.contigs(file=data/WAVES_data/mothur/WAVES.file, outputdir=data/WAVES_data/mothur/, processors=2)
```

Check out sequences, see how they look

```
summary.seqs(fasta=current)
```

  Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	239	239	0	3	1
2.5%-tile:	1	253	253	0	4	9218
25%-tile:	1	253	253	0	4	92176
Median: 	1	253	253	0	4	184352
75%-tile:	1	253	253	0	5	276528
97.5%-tile:	1	424	424	4	9	359486
Maximum:	1	500	500	50	234	368703
Mean:	1	259.007	259.007	0.42713	4.73255
Number of Seqs:	368703

So, looks like we have a few bad ones. Need to do the filtering that Pat showed us. The long polymer is especially bad.

```
screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=80)

```
Output
##########
mothur > summary.seqs(fasta=current)
Using data/WAVES_data/mothur/seq.pairs.trim.contigs.good.fasta as input file for the fasta parameter.

Using 2 processors.

		Start	End	NBases	Ambigs	Polymer	NumSeqs
Minimum:	1	250	250	0	3	1
2.5%-tile:	1	253	253	0	4	7756
25%-tile:	1	253	253	0	4	77559
Median: 	1	253	253	0	4	155117
75%-tile:	1	253	253	0	5	232675
97.5%-tile:	1	254	254	0	6	302477
Maximum:	1	273	273	0	11	310232
Mean:	1	253.024	253.024	0	4.50617
No. of Seqs:	310232

##########

Still have a homopolymer of 11 for some reason. But, not looking too bad. Let's just move ahead.

Grab unique sequences.

```
unique.seqs(fasta=current)
```


Prep for next time

[WARNING]: Cannot determine number of physical pages

Current RAM usage: 0.0804634 Gigabytes. Total Ram: 0 Gigabytes.

Current files saved by mothur:
fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.fasta
group=data/WAVES_data/mothur/WAVES.contigs.good.groups
name=data/WAVES_data/mothur/WAVES.trim.contigs.good.names
qfile=data/WAVES_data/mothur/WAVES.trim.contigs.qual
processors=2
summary=data/WAVES_data/mothur/WAVES.trim.contigs.good.summary

Current input directory saved by mothur: data/WAVES_data/mothur/

Current output directory saved by mothur: data/WAVES_data/mothur/

Current default directory saved by mothur: /home/ubuntu/mothur/

Current working directory: /home/ubuntu/

Output File Names:
data/WAVES_data/mothur/current_files.summary



#### 07-17-16

Follow steps 5 and 6 and seven to connect to ec2 (with different address) with the volume

Continue on with mothur pipeline

Change "current" settings

```
set.current(fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.fasta, group=data/WAVES_data/mothur/WAVES.contigs.good.groups, name=data/WAVES_data/mothur/WAVES.trim.contigs.good.names, qfile=data/WAVES_data/mothur/WAVES.trim.contigs.qual, processors=2, summary=data/WAVES_data/mothur/WAVES.trim.contigs.good.summary, inputdir=data/WAVES_data/mothur/, outputdir=data/WAVES_data/mothur/)
```

Just made the unique seqs files, so now we need to know how many of those unique seqs there are. Remember, this is all to simplify the amount of computation that goes into the final alignment.

```
count.seqs(name=WAVES.trim.contigs.good.names, group=current)
```

Now we have a file (WAVES.trim.contigs.good.count_table) that has the number of sequences that fit each unique contig.

Now we need to align this to the known region of the 16S sequence. This can be done with other types of sequences, if we were interested in amplicon sequencing (which I will be!) Start by excising the 16S gene from the silva database.

#### 07-18-16 Continue mothur analysis

Pull out 16S sequences from Silva database.

```
pcr.seqs(fasta=data/references/silva.bacteria/silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=2)
```

Rename it to a more informative name

```
system(mv data/WAVES_data/mothur/silva.bacteria.pcr.fasta data/WAVES_data/mothur/silva.v4.reference.fasta)
```

Check the file that we just made

```
summary.seqs(fasta=silva.v4.reference.fasta)
```

Now that we have a reference alignment, we can align our sequences to the silva v4 reference.

```
align.seqs(fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.fasta, reference=data/WAVES_data/mothur/silva.v4.reference.fasta)
```

Needed to change current for count

```
set.current(count=WAVES.trim.contigs.good.count_table)
```

Get a summary table. The results show us that not all of our sequences align with the regions that we would expect it to (which would be 1968 to 11550). So now we need to trim that.

```
summary.seqs(fasta=current, count=current)

#Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	3	1	0	1	1
#2.5%-tile:	1968	11550	253	0	4	7756
#25%-tile:	1968	11550	253	0	4	77559
#Median: 	1968	11550	253	0	4	155117
#75%-tile:	1968	11550	253	0	5	232675
#97.5%-tile:	1968	11550	254	0	6	302477
#Maximum:	13425	13425	272	0	11	310232
#Mean:	1968.97	11548.9	252.983	0	4.50565
## of unique seqs:	36936
#total # of seqs:	310232
##############################
```

Trimming the alignments. Make sure to trim both the unique OTU fasta file and the count table
```
screen.seqs(fasta=current, count=current, summary=current, start=1968, end=11550, maxhomop=8)
```

Double-check the sequences again to make sure it trimmed them properly.

```
summary.seqs(fasta=current, count=current)

#Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	11550	250	0	3	1
#2.5%-tile:	1968	11550	253	0	4	7725
#25%-tile:	1968	11550	253	0	4	77248
#Median: 	1968	11550	253	0	4	154496
#75%-tile:	1968	11550	253	0	5	231744
#97.5%-tile:	1968	11550	254	0	6	301267
#Maximum:	1968	13406	272	0	8	308991
#Mean:	1967.97	11550	253.025	0	4.50503
## of unique seqs:	36529
#total # of seqs:	308991
#
#Output File Names:
#data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.summary
#
#It took 10 secs to summarize 308991 sequences.
```

Filter it so that we only have the ~250bp sequence. The "trump=." argument removes columns that have a . in that column for ANY sequence.

```
 filter.seqs(fasta=current, vertical=T, trump=.)
 ```

Now that we've filtered the sequences, there will likely be a few more sequences that are similar, and can combine them using unique.seqs

```
unique.seqs(fasta=current, count=current)
```

Check it! Looks like we reduced the number of unique samples by 12 after the latest round of filtering.

```
summary.seqs(fasta=current)

#		Start	End	NBases	Ambigs	Polymer	NumSeqs
#Minimum:	1	482	233	0	3	1
#2.5%-tile:	1	483	252	0	4	913
#25%-tile:	1	483	253	0	4	9130
#Median: 	1	483	253	0	4	18259
#75%-tile:	1	483	253	0	5	27388
#97.5%-tile:	1	483	254	0	6	35605
#Maximum:	1	483	265	0	8	36517
#Mean:	1	483	253.005	0	4.48742
## of Seqs:	36517
#
#Output File Names:
#data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.#unique.summary
```

mothur > get.current()
[WARNING]: Cannot determine number of physical pages

Current RAM usage: 0.226707 Gigabytes. Total Ram: 0 Gigabytes.

Current files saved by mothur:
fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.fasta
group=data/WAVES_data/mothur/WAVES.contigs.good.groups
name=data/WAVES_data/mothur/WAVES.trim.contigs.good.names
qfile=data/WAVES_data/mothur/WAVES.trim.contigs.qual
count=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.count_table
processors=2
summary=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.summary

Current input directory saved by mothur: data/WAVES_data/mothur/

Current output directory saved by mothur: data/WAVES_data/mothur/

Current default directory saved by mothur: /home/ubuntu/mothur/

Current working directory: /home/ubuntu/

Output File Names:
data/WAVES_data/mothur/current_files.summary


#### 07-25-16 Continue mothur analysis

Load ec2 instance with previous parameters, mount volume, load up mothur and change current to above settings

```
set.current(fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.fasta, group=data/WAVES_data/mothur/WAVES.contigs.good.groups, name=data/WAVES_data/mothur/WAVES.trim.contigs.good.names, qfile=data/WAVES_data/mothur/WAVES.trim.contigs.qual, count=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.count_table, processors=2, summary=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.summary, inputdir=data/WAVES_data/mothur/, outputdir=data/WAVES_data/mothur/)
```

Pre-cluster the data. This is to denoise the data set. It works by combining sequences that have x number of amino acid differences. x is set by the command "diffs=". Pat suggested setting this to 1 per 100bp you have in your sequence. We have 250bp sequences, so we'll use diffs=2. I'd rather err on the side of pre-clustering too little than too much. Because I know oh so much about how this works.

```
pre.cluster(fasta=current, count=current, diffs=2)

#Total number of sequences before pre.cluster was 36517.
#pre.cluster removed 28934 sequences.
#
#It took 14 secs to cluster 36517 sequences.
#It took 16 secs to run pre.cluster.
#
#Output File Names:
#data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.fasta
#data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.count_table
#data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.WAVES.map
```

Note: We removed a shit-ton of sequences. But, in the MiSeq SOP with Pat, we reduced the number of sequences from 16349 to 5657. So, probably nothing to worry about. Might want to look into this at some point though, and what it means.

Now to look for chimeras. We set dereplicate=T to prevent the script from removing sequences from one sample that were flagged as chimeras in another sample. Pat suggested that the default of dereplicate=F is too aggressive.

```
chimera.uchime(fasta=current, count=current, dereplicate=t)

04:55 7.2Mb  100.0% 979/7582 chimeras found (12.9%)

It took 295 secs to check 7583 sequences from group WAVES.

Output File Names:
data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table
data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.chimeras
data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos
```

This only removes chimeras from the count table, but it does store the sequences that it removed so that we can remove them from the fasta file. Here's how we do that:

```
remove.seqs(fasta=current, accnos=WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos)

# Removed 979 sequences from your fasta file.
#
#Output File Names:
#data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
```

NOTE: Could have used accnos=current here, that is saved as a current. Will need to do this if I make a workflow of it.

Okay. So there were 979 chimeras removed from the fasta files. The final number of sequences was 307085, which is down from 308991. These means that 12.9% of the unique sequences were identified as chimeras, but only 0.42% of all sequences were identified as chimeras. This seems low, in the MiSeq SOP dataset they remove ~7% of the sequences as chimeras. I don't know how to identify ours as wrong, though, so we'll go with it for now. This is something to look back into though.

Now that chimeras are gone, we should remove unwanted lineages. Will do this using the RDP training set. Needed to download them and upload them to our Volume, since we're not in the mothur AMI as "mothur". Pretty simple to download though. Will quit mothur, and upload them.

current for when we go back:
__________________________________________________

Current files saved by mothur:
accnos=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos
fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta
group=data/WAVES_data/mothur/WAVES.contigs.good.groups
name=data/WAVES_data/mothur/WAVES.trim.contigs.good.names
qfile=data/WAVES_data/mothur/WAVES.trim.contigs.qual
count=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table
processors=2
summary=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.summary

Current input directory saved by mothur: data/WAVES_data/mothur/

Current output directory saved by mothur: data/WAVES_data/mothur/

Current default directory saved by mothur: /home/ubuntu/mothur/

Current working directory: /home/ubuntu/

Output File Names:
data/WAVES_data/mothur/current_files.summary
__________________________________________________


Goddammit, I'm an idiot. All the references, data, code, etc, is still in the instance if you don't log in as mothur, it's just up a directory. No need to upload the data training set, I don't think. Let's see if we can still access them.

/home/mothur/data/references/
```
set.current(accnos=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos, fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, group=data/WAVES_data/mothur/WAVES.contigs.good.groups, name=data/WAVES_data/mothur/WAVES.trim.contigs.good.names, qfile=data/WAVES_data/mothur/WAVES.trim.contigs.qual, count=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, processors=2, summary=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.summary, inputdir=data/WAVES_data/mothur/, outputdir=data/WAVES_data/mothur)
```

Okay. To remove unwanted lineages:

```
classify.seqs(fasta=current, count=current, reference=/home/mothur/data/references/trainset14_032015.pds.fasta, taxonomy=/home/mothur/data/references/trainset14_032015.pds.tax, cutoff=80)
```

Didn't work. Try copying the references into my own reference folder

```
cp ../mothur/data/references/trainset14_032015.pds.tax data/references/
cp ../mothur/data/references/trainset14_032015.pds.fasta data/references/

set.current(accnos=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos, fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, group=data/WAVES_data/mothur/WAVES.contigs.good.groups, name=data/WAVES_data/mothur/WAVES.trim.contigs.good.names, qfile=data/WAVES_data/mothur/WAVES.trim.contigs.qual, count=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, processors=2, summary=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.summary, inputdir=data/WAVES_data/mothur/, outputdir=data/WAVES_data/mothur)
```

Now try to classify them using the database

```
classify.seqs(fasta=current, count=current, reference=data/references/trainset14_032015.pds.fasta, taxonomy=data/references/trainset14_032015.pds.tax, cutoff=80)
```

Alright, it's already working. Cool cool cool. Might take awhile. Here's the output:

```
#It took 145 secs to classify 6604 sequences.
#
#It took 0 secs to create the summary file for 6604 sequences.
#
#Output File Names:
#data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy
#data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary
```

Now remove mitochondria, chloroplast, eukaryote, and unknown sequences
```
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Eukaryota)
```

Whoa. Removed a lot of sequences (according to summary.seq). These numbers are down from 6604 and 307085

```
# of unique seqs:	6011
total # of seqs:	244802
```

That's 9% of unique sequences removed, and 20.3% of total sequences. Can't be right, right? Keep going now, get to our classifications, check it later if it seems weird. Which it seems like it might.

Now, we cluster to OTUs and phylotypes. Use the cluster.split functions, which splits the data into groups based on phylogenetic classification, then clusters them from there into OTUs. It's faster than more traditional clustering, and shouldn't change much in the final result.

```
cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.15)
```

Whoa, very fast.

How many sequences are from each OTU in each group?

```
make.shared(list=current, count=current, label=0.03)
```

This generates a shared table, which represents the number of times an OTU is seen in each table.

We want to know the taxonomic classification of each OTU, so we generate a consensus classification for each OTU.

```
classify.otu(list=current, count=current, taxonomy=current, label=0.03)
```

Stop here, as far as data-manipulation goes. We could also group the sequences by phylotype or make a phylogenetic tree, but we'll head straight into the analysis. I think we may want to re-run this all with green genes anyways. If I do that, I'll check out Robin's workflow, and try running this data through that.

Check out the number of counts from each sample

```
count.groups(shared=current)
```

Uh-oh. This appears to be returning one count, it's using WAVES as a grouping. This needs to be done for each individual sample I think. This goes way back to the beginning.

Here's current, just in case.

_____________________________________________
Current files saved by mothur:
accnos=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos
column=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta.94.dist
fasta=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta
group=data/WAVES_data/mothur/WAVES.contigs.good.groups
list=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list
name=data/WAVES_data/mothur/WAVES.trim.contigs.good.names
qfile=data/WAVES_data/mothur/WAVES.trim.contigs.qual
shared=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared
taxonomy=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
count=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table
processors=2
summary=data/WAVES_data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.pick.summary

Current input directory saved by mothur: data/WAVES_data/mothur/

Current output directory saved by mothur: data/WAVES_data/mothur/

Current default directory saved by mothur: /home/ubuntu/mothur/

Current working directory: /home/ubuntu/

_____________________________________________


Fuck! It sees the underscore in the WAVES_data directory and is using that as the group name! Might need to stop naming my directories with underscores?

Damn. Even changing into the directory doesn't help. I need to change the directory names. FUCK. Changed underscores in the path to periods. Changed mothur directory in WAVEs.data to mothur.old, created new mothur directory to work with. Here's the code.

Separate code for making the 16S alignment file from Silva database. Basically copied from the first iteration of this attempt. (I actually just copied the old file into the new mothur file, but wanted the code that generated it shown here).

```
# Pull out the correct sequence alignment from the Silva database to generate a file against which we can align our 16S genes
pcr.seqs(fasta=data/references/silva.bacteria/silva.bacteria.fasta, outputdir=data/WAVES.data/mothur/, start=11894, end=25319, keepdots=F, processors=2)

# Rename it to a more informative name
system(mv data/WAVES.data/mothur/silva.bacteria.pcr.fasta data/WAVES.data/mothur/silva.v4.reference.fasta)
```

###### Mothur workflow to clean up the sequences.

```
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
```

###### Sanity checks!

Check the alignment

```
summary.seqs(fasta=data/WAVES.data/mothur/WAVES.trim.contigs.good.unique.align, count=data/WAVES.data/mothur/WAVES.trim.contigs.good.count_table)
```

How much did pre-clustering remove?

```
summary.seqs(fasta=data/WAVES.data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.fasta, count=data/WAVES.data/mothur/WAVES.trim.contigs.good.unique.good.filter.count_table)

# vs.

summary.seqs(fasta=data/WAVES.data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=data/WAVES.data/mothur/WAVES.trim.contigs.good.unique.good.filter.unique.precluster.count_table)
```
