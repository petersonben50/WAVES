# Always get a clean start
rm(list = ls())

# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####

# Load up needed packages
##########################
#####

# Read in necessary files and alter them as necessary
##########################
# Read in the relatived OTU table.
WAVES.REL <- read.table("dataEdited/WAVES.final.abund",
                        # OTU numbers are the heading names
                        header = TRUE,
                        # The file is tab-delimited
                        sep = "\t",
                        # Use the site names as row names
                        row.names = 1)
# Convert WAVES.REL from a dataframe to a matrix
WAVES.REL <- as.matrix(WAVES.REL)
# Transpose matrix
WAVES.REL <- t(WAVES.REL)

# Generate a matrix with the taxonomic classifications for each OTU.
# Read in OTU taxonomy table
WAVES.tax <- read.csv2("dataEdited/WAVES.taxonomy",
                       stringsAsFactors = FALSE,
                       header = FALSE)
# Sort OTU tax table by OTU #
WAVES.tax <- WAVES.tax[order(WAVES.tax$V1), ]
# Make OTU names the row names
otu.names <- WAVES.tax[, 1]
rownames(WAVES.tax) <- otu.names
WAVES.tax <- WAVES.tax[, -1]
# Name the columns for the phylogenetic classification
phylo.names <- c("Kingdom",
                 "Phylum",
                 "Class",
                 "Order",
                 "Family",
                 "Genus",
                 "Species")
colnames(WAVES.tax) <- phylo.names
WAVES.tax <- WAVES.tax[, phylo.names]
# Remove bootstrap values
for (i in seq_len(ncol(WAVES.tax))) {
  WAVES.tax[, i] <- sapply(WAVES.tax[, i], 
                                function(x) {
                                  sapply(strsplit(x,
                                                  "\\("),
                                         "[", 1)
                                })
}
# Remove the tax. classification level 
for (i in 1:length(WAVES.tax)) {
  for (j in grep("[k,p,f,c,o,g]*__", WAVES.tax[[i]])) {
    WAVES.tax[[i]][j] <- strsplit(WAVES.tax[[i]][j],
                                  "[k,p,f,c,o,g]*__")[[1]][2]
  }
}


# Read in list of methylators
methylators <- read.csv("../metagHgcA/references/hgcA_compiled_list.csv",
                        stringsAsFactors = FALSE)
methylators <- methylators[, 1:3]
meth.genus <- sapply(methylators$strain, 
       function(x) {
         strsplit(x,
                  " ")[[1]][1]
})
meth.species <- sapply(methylators$strain, 
                     function(x) {
                       strsplit(x,
                                " ")[[1]][2]
                     })
#####


# Look for taxonomic groups that might be methylators
# Not a very good way to look for them, I admit.
##########################
# Let's start easy. Are there Geobacter? Where are they?
WAVES.REL[grep("Geobacter",
               WAVES.tax$Genus,
               ignore.case = TRUE), ]
# Only at 19m on Day B

# How about desulfos?
WAVES.REL[grep("Desulfovibrio",
               WAVES.tax$Genus,
               ignore.case = TRUE), ]
# Again, only see them very deep.

# How about "all" (code is messy) methylating genera?
lapply(meth.species, function(x) {
  grep(x,
       WAVES.tax[, 7],
       ignore.case = TRUE)
})