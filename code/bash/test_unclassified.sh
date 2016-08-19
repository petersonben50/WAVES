grep 'Bacteria_unclassified' dataEdited/WAVES.taxonomy > testing/unclassified.txt

tr ';k' '
' < testing/unclassified.txt | grep "Otu" > testing/list.unclassified

grep -f testing/list.unclassified /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.below.98.all | wc -l


# For non-representative sequences

# Move all sequences unclassified at the phylum level into a file.
# This is for use with classification of all OTUs, not just the
# representative sequences.
#grep -nr "M0.*;.*;k__Bacteria_unclassified.*;.*;.*;.*;.*;.*" dataEdited/WAVES.all.taxonomy > testing/all.unclassified.txt

# Make a file that contains a list of all the sequences that are
# unclassified at the phylum level.
#tr ';k' '
#' < testing/all.unclassified.txt | grep "M0" > testing/test.txt
# Remove first part of file.
#tr ':' '
#' < testing/test.txt | grep "M0" > testing/list.unclassified.txt

# How many sequences from ids.below.98.all are in that group of unclassified sequences?
#grep -f testing/list.unclassified.txt /Volumes/BDPLabHD/data/WAVES/dataEdited/16STaxAss/ids.below.98.all | wc -l
