# This code generates a heatmap that includes the top 
# twenty orders from each sample.


# Always get a clean start
rm(list = ls())


# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####


# Load up needed packages
##########################
library("ggplot2")
library("gplots")
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


# Generate a heatmap that includes the top 10 families
# in each sample. Obviously many of these will be redundant.
##########################
# Combine by family
WAVES.REL.family <- aggregate(cbind(A1002, A1009, A1010, A1011, A1012, A1013, A1216, B1002, B1012, B1014, B1015, B1016, B1019) ~ WAVES.tax$Family,
                               data = WAVES.REL,
                               FUN = sum)
# Make the family names the row names
rownames(WAVES.REL.family) <- WAVES.REL.family[[1]]
WAVES.REL.family <- WAVES.REL.family[-1]
# Make this a data matrix
WAVES.REL.family <- as.matrix(WAVES.REL.family)
# Remove the "f__" from the name
for (i in grep("f__", rownames(WAVES.REL.family))) {
  rownames(WAVES.REL.family)[i] <- strsplit(rownames(WAVES.REL.family)[i],
                                "f__")[[1]][2]
}

# Pull out the top 10 orders in each grouping.
WAVES.families <- matrix(0,
                         nrow = 10,
                         ncol = 13)
for (i in seq_len(ncol(WAVES.REL.family))) {
  WAVES.families[, i] <- row.names(WAVES.REL.family)[head(sort(WAVES.REL.family[, i],
                                                  decreasing = TRUE,
                                                  index.return = TRUE)$ix,
                                             10)]
  }
# Find which ones are unique
WAVES.families <- unique(as.vector(WAVES.families))
# Remove unclassified samples
WAVES.families <- WAVES.families[-grep("unclassified", WAVES.families)]
# Remove NA values
WAVES.families <- WAVES.families[!is.na(WAVES.families)]

# Make the heatmap
heatmap.2(apply(WAVES.REL.family[WAVES.families,],
                2,
                sqrt),
                 #col=hc(20),
                 #scale="col",
                 key=TRUE,
                 symkey=FALSE,
                 trace="none",
                 density.info="none",
                 dendrogram="both",
                 #labCol = group.for.comparison,
                 #colCol=c(map$Class_Color),
                 margins=c(5,13), srtCol=90)
#####