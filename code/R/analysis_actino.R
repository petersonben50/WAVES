# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####

# Load up needed packages
##########################
library("ggplot2")
library("RColorBrewer")
library("gplots")
#####

# Load in data and transform it
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

# Generate taxonomy data.fram
WAVES.tax.list <- matrix(data = 0,
                         nrow = 3387,
                         ncol = 6)
# Remove bootstrap values
for (i in seq_len(ncol(WAVES.tax))) {
  WAVES.tax.list[, i] <- sapply(WAVES.tax[, i], 
                                function(x) {
                                  sapply(strsplit(x,
                                                  "\\("),
                                         "[", 1)
                                })
}
# Add in phylogenetic levels for names
colnames(WAVES.tax.list) <- names(WAVES.tax)
# Pull out just the actinobacteria at the phylum level
WAVES.tax.actino <- WAVES.tax.list[WAVES.tax.list[, "Phylum"] == "p__Actinobacteria", ]
WAVES.REL.actino <- WAVES.REL[WAVES.tax.list[, "Phylum"] == "p__Actinobacteria", ]
#####

# Family analysis
##########################
# Combine the bacteria by genus
WAVES.REL.actino.genus <- aggregate(cbind(A1002, A1009, A1010, A1011, A1012, A1013, A1216, B1002, B1012, B1014, B1015, B1016, B1019) ~ WAVES.tax.actino[, "Genus"],
                               data = WAVES.REL.actino,
                               FUN = sum)
# Name the columns by genus
rownames(WAVES.REL.actino.genus) <- unique(WAVES.tax.actino[, "Genus"])
WAVES.REL.actino.genus <- WAVES.REL.actino.genus[, -1]
# Make the heatmap
fig5B<-heatmap.2(as.matrix(WAVES.REL.actino.genus),
                 #col=hc(20),
                 scale="none",
                 key=TRUE,
                 symkey=FALSE,
                 trace="none",
                 density.info="none",
                 dendrogram="both",
                 #labCol = group.for.comparison,
                 #colCol=c(map$Class_Color),
                 margins=c(5,13), srtCol=90)
#####