# This code generates a rank abundance curve for
# each of our samples


# Always get a clean start
##########################
rm(list = ls())
#####

# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####


# Load up needed packages
##########################
library("reshape2")
library("vegan")
library("BiodiversityR")
#####


# Read in necessary files
##########################
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
WAVES.tax <- WAVES.tax[, phylo.names]
# Remove the bootstrap values from the phylogenetic classification
for (i in seq_len(ncol(WAVES.tax))) {
  WAVES.tax[, i] <- sapply(WAVES.tax[, i], 
                           function(x) {
                             sapply(strsplit(x,
                                             "\\("),
                                    "[", 1)
                           })
}

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

# Establish site names
site.names <- colnames(WAVES.REL)
#####


##########################
# Results
##########################


# Make a rank abundance curve of each sample. 
# Run at phylum level.
##########################

# First, pull out phylum level abundances. 
# Borrowed from code in analysis_general.R (from this project)
# Slight modifications

    # ***Requires "reshape2"***
    # First, make a vector with the phyla names that correspond
    # to the OTU names
    # This splits the phylum classification at the __, then splits 
    # the second element of that first list at the (. Finally,
    # the first part of that last list is returned. This gives us
    # the vector we need.
    WAVES.phyla <- sapply(strsplit(WAVES.tax$Phylum, 
                                   "__"), 
                          function(x) {
                            sapply(strsplit(x[2],
                                            "\\("),
                                   "[", 1)
                          })
    # Make a new abundance table with the phylum names for the rows
    WAVES.REL <- as.data.frame(WAVES.REL)
    WAVES.REL.phyla <- cbind(WAVES.phyla, WAVES.REL)
    # Sum the rows that are of the same phyla
    WAVES.REL.phyla <- aggregate(cbind(A1002, A1009, A1010, A1011, A1012, A1013, A1216, B1002, B1012, B1014, B1015, B1016, B1019) ~ WAVES.phyla,
                                 data = WAVES.REL.phyla,
                                 FUN = sum)
    # Melt the data
    melt.WAVES.REL.phyla <- melt(WAVES.REL.phyla, 
                                 id.vars=c("WAVES.phyla"))
    # If the phyla never accounts for more than 0.5% of the sample, remove them
    # Convert them to 0 if they are less than 0.005
    melt.WAVES.REL.phyla[melt.WAVES.REL.phyla$value <= 0.005, "value"] <- 0
    # Recast it so that the columns are the samples, the rows are the phyla
    # Cell values are the relativized abundance
    WAVES.REL.recast <- dcast(melt.WAVES.REL.phyla, 
                              formula = WAVES.phyla ~ variable)
    # Remove phyla with a row sum of 0.
    WAVES.REL.recast <- WAVES.REL.recast[rowSums(WAVES.REL.recast[, -1]) != 0, ]
# The rank abundance function requires phyla in columns, samples in rows.
# First need to set phyla as row names
row.names(WAVES.REL.recast) <- WAVES.REL.recast[, 1]
# Transpose the matrix without the names column
WAVES.REL.phyla.ra <- t(WAVES.REL.recast[, -1])

# Rank abundance curve
# Just make a fucking barplot rather than trying to mess with BiodiversityR
barplot(sort(WAVES.REL.phyla.ra[1, ],
             decreasing = TRUE),
        axisnames = FALSE)
axis(side = 1,
     at = seq(from = 0.6,
              to = 14.5,
              by = 1.22),
     names(sort(WAVES.REL.phyla.ra[1, ],
           decreasing = TRUE)
           ),
     las = 2,
     padj = 0.5
     )
points(x = 0.6,
       y = 0.05)      