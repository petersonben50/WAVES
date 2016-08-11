#### Classification of WAVES sequences using Robin's 16STaxAss workflow

The sequences and the reference databases are stored on BDPLabHD. The sample FASTA files I am working with have already been aligned and quality-filtered using mothur. The output files

###### First, change into WAVES project directory

```
cd Documents/gradSchool/research/mercury/WAVES/
```

###### Generate a GreenGenes database that is compliant with the needs for this workflow.

They need to match the freshwater database formatting. Let's copy a new file so we don't lose the original formatting.

```
cp /Volumes/BDPLabHD/data/referenceDB/gg_13_5_taxonomy.txt /Volumes/BDPLabHD/data/referenceDB/gg_13_5_taxonomy_FWedit.txt
```

Remove spaces in the copied Green Genes database, store as NoSpaces in working directory
```
sed 's/ //g' </Volumes/BDPLabHD/data/referenceDB/gg_13_5_taxonomy_FWedit.txt >NoSpaces
```

Add semicolon to end of each line of the new NoSpaces file, save as EndLineSemicolons
```
sed 's/$/;/' <NoSpaces >EndLineSemicolons
```

Return edited file to the referenceDB directory on the hard drive
```
mv EndLineSemicolons /Volumes/BDPLabHD/data/referenceDB/gg_13_5_taxonomy_FWedit.txt
```

Remove the NoSpaces file
```
rm NoSpaces
```

The fasta file is okay though, I think.

###### Blast our sequences against the FW database.

I want to use the list of unique, filtered sequences that I generated using mothur's quality control measures. I should be able to use this to classify my OTUs.

Going to run this on the filterTest data first, see where I come out.
