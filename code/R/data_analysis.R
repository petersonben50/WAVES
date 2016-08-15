# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####

# Load up needed packages
##########################
library("vegan")
#####

# Read in necessary files
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
#####


# Generate a PCoA graph from the relativized data.
##########################
# Generate a distance matrix between the samples using the weighted data.
WAVES.REL.dist <- vegdist(WAVES.REL,
                          method = "bray",
                          )