# Always get a clean start
rm(list = ls())


# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####


# Read in necessary files and alter them as needed.
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

# Generate a matrix that contains the taxonomic information for each OTU
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
                 "Genus")
colnames(WAVES.tax) <- phylo.names
# Keep the columns only to the Genus/Lineage level.
# Don't want/need Species or Tribe
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
#####

##########################
# Analysis
##########################


# How abundant are Chlorobi in each sample?
##########################
# Aggregate relativized abundance by phyla
WAVES.REL.phyla <- aggregate(cbind(A1002, A1009, A1010, A1011, A1012, A1013, A1216, B1002, B1012, B1014, B1015, B1016, B1019) ~ WAVES.tax[, "Phylum"],
                             data = WAVES.REL,
                             FUN = sum)
# Rownames as phyla
rownames(WAVES.REL.phyla) <- WAVES.REL.phyla[, "WAVES.tax[, \"Phylum\"]"]
# Remove column with names
WAVES.REL.phyla <- WAVES.REL.phyla[, -1]
# Make into matrix
WAVES.REL.phyla <- as.matrix(WAVES.REL.phyla)
# Pull out Chlorobi
WAVES.REL.phyla["Chlorobi", ]
# Meh. Tough to tell since this is all relativized data anyways.
#####


# What families of Chlorobi are present?
##########################
chlorobi.fam.abund <- cbind(WAVES.tax[which(WAVES.tax$Phylum == "Chlorobi"), "Family"],
                            WAVES.REL[which(WAVES.tax$Phylum == "Chlorobi"), ])
#####