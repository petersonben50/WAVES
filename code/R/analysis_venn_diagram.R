# Venn Diagram of the OTUs
##########################
# I'm going to create a Venn diagram from our sample. I will 
# cluster the samples first and combine those that clustered together
# so that I can create a Venn diagram of 4-5 different "groups".
# Maybe I should Venn diagram the samples in the epilimnion to 
# see if they really are well mixed.


# Start with a clean slate.
rm(list = ls())

# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####

# Load up needed packages
##########################
library("vegan")
library("gplots")
#library("reshape2")
#library("ggplot2")
#library("scales")
#library("RColorBrewer")
#####

# Read in relativized abundance table and alter it.
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
site.names <- colnames(WAVES.REL)


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
# Remove bootstrap values from names
for (i in seq_len(ncol(WAVES.tax))) {
  WAVES.tax[, i] <- sapply(WAVES.tax[, i], 
                           function(x) {
                             sapply(strsplit(x,
                                             "\\("),
                                    "[", 1)
                           })
}
#####

# Read in the environmental data
##########################
# Read in the .csv file with the environmental data
# Transpose it so that the samples are the rows.
WAVES.meta <- as.data.frame(t(read.csv("dataRaw/2,30 Sept 2009 ME canoco.csv",
                                       header = TRUE,
                                       row.names = 1,
                                       stringsAsFactors = FALSE)))
# Keep only the values for the samples that we sequenced.
WAVES.meta <- WAVES.meta[colnames(WAVES.REL), ]
# Store the site names in a vector for later.
site.names <- row.names(WAVES.meta)
# Remove spaces in the column headers.
colnames(WAVES.meta) <- gsub(" ", 
                             ".", 
                             colnames(WAVES.meta))
# Pull out environmental data that might be relevant
WAVES.meta <- WAVES.meta[, c("Day",
                             "Depth",
                             "Temp",
                             "DO.(mg/L)",
                             "pH")]
# Convert all columns to numerics
WAVES.meta[, -1] <- apply(WAVES.meta[, -1],
                          2,
                          as.numeric)
# Add in numeric dates for time points.
WAVES.meta[, "Day"] <- gsub("A", 
                            "2009/09/02",
                            WAVES.meta[, "Day"])
WAVES.meta[, "Day"] <- gsub("B", 
                            "2009/09/30",
                            WAVES.meta[, "Day"])
# Need to add in row names again.
row.names(WAVES.meta) <- site.names
#####

# Generate a distance matrix between the samples using the weighted data.
##########################
# Use vegan to generate a Bray-Curtis distance matrix
WAVES.REL.dist <- vegdist(t(WAVES.REL),
                          method = "bray")
# Generate the PCoA from the distance matrix
WAVES.REL.pcoa <- cmdscale(WAVES.REL.dist)
# Set up the PCoA plot but, don't include the data yet.
plot(WAVES.REL.pcoa[, 1],
     WAVES.REL.pcoa[, 2],
     cex = 0,
     main = "WAVES data PCoA")
# Add data to the plot.
# Points are names on the plot
text(WAVES.REL.pcoa[, 1],
     WAVES.REL.pcoa[, 2],
     attr(WAVES.REL.pcoa, "dimnames")[[1]],
     cex = 0.6)
# This makes it look like I could use 5 "groups" for the
# venn diagram. I worry about doing this in that it wouldn't 
# be even since the 2 epilimnetic groups would have more
# samples.
#####

# Venn diagram of the 6 samples in the epilimnion of 
# 9-2-09 that clustered together
# Interesting that A-13 is above the thermocline, but below the oxycline. 
##########################
A.epi <- c("A1002", "A1009", "A1010",
           "A1011", "A1012", "A1013")
# I think we can only do 4. Let's cluster these 6
# sequences and see if we can get it down to 4 groups.
WAVES.REL.A.epi.dist <- vegdist(t(WAVES.REL[, A.epi]),
                                method = "bray")
# Generate the PCoA from the distance matrix
WAVES.REL.A.epi.pcoa <- cmdscale(WAVES.REL.A.epi.dist)
# Set up the PCoA plot but, don't include the data yet.
plot(WAVES.REL.A.epi.pcoa[, 1],
     WAVES.REL.A.epi.pcoa[, 2],
     cex = 0,
     main = "WAVES data PCoA")
# Add data to the plot.
# Points are names on the plot
text(WAVES.REL.A.epi.pcoa[, 1],
     WAVES.REL.A.epi.pcoa[, 2],
     attr(WAVES.REL.A.epi.pcoa, "dimnames")[[1]],
     cex = 0.6)
# Not that surprising. The samples get further and further
# apart the lower they go in the epilimnion. The top most
# three cluster.
#####

# Let's try a Venn diagram of A-2, A-16, B-2, and B-16
##########################
# Now make it a presence/absence table.
WAVES.PA <- (WAVES.REL > 0)*1
# There are the four samples of interest.
sites <- c("A1002", "B1002", "B1016", "A1216")
# Need a presence/absence table for the venn function
# of gplots
WAVES.PA.venn <- WAVES.PA[, sites]
venn(WAVES.PA.venn)
#####

# Venn diagram of B samples
##########################
B.epi <- c("B1002", "B1012", "B1014",
           "B1015", "B1016", "B1019")
# I think we can only do 4. Let's cluster these 6
# sequences and see if we can get it down to 4 groups.
WAVES.REL.B.epi.dist <- vegdist(t(WAVES.REL[, B.epi]),
                                method = "bray")
# Generate the PCoA from the distance matrix
WAVES.REL.B.epi.pcoa <- cmdscale(WAVES.REL.B.epi.dist)
# Set up the PCoA plot but, don't include the data yet.
plot(WAVES.REL.B.epi.pcoa[, 1],
     WAVES.REL.B.epi.pcoa[, 2],
     cex = 0,
     main = "WAVES data PCoA")
# Add data to the plot.
# Points are names on the plot
text(WAVES.REL.B.epi.pcoa[, 1],
     WAVES.REL.B.epi.pcoa[, 2],
     attr(WAVES.REL.B.epi.pcoa, "dimnames")[[1]],
     cex = 0.6)
# Not that surprising. The samples get further and further
# apart the lower they go in the epilimnion. The top most
# three cluster.
B.epi <- c("B1002", "B1012", "B1014", "B1015")
venn((WAVES.REL[, B.epi] >0)*1)
#####

