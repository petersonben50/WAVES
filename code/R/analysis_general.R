# Always get a clean start
rm(list = ls())

# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####

# Load up needed packages
##########################
library("vegan")
library("phyloseq")
library("reshape")
library("reshape2")
library("ggplot2")
library("scales")
library("RColorBrewer")
library("gplots")
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

# Read in buoy profile
temp.buoy.02 <- t(read.csv("dataRaw/20090902_buoy.csv",
                      stringsAsFactors = FALSE,
                      header = FALSE,
                      nrows = 1))
temp.buoy.30 <- t(read.csv("dataRaw/20090930_buoy.csv",
                        stringsAsFactors = FALSE,
                        header = FALSE,
                        nrows = 1))
temp.buoy <- cbind(temp.buoy.02,
                   temp.buoy.30)
# Keep the temp measurements
temp.buoy <- temp.buoy[11:33, ]
# Add column with depths
temp.buoy <- cbind(temp.buoy,
                   0:22)
# Name the columns
colnames(temp.buoy) <- c("2009/09/02",
                         "2009/09/30",
                         "Depth")
#####

# Overview of Environmental Data
##########################
plotstuff <- function(yearmonthday) {
  # Set up plot parameters
  par(mar = c(6, 3, 3, 2),
      mgp = c(3, 0.5, 0))
  
  # Set up plot, put temperature data points on it.
  plot(temp.buoy[, yearmonthday],
       -temp.buoy[, "Depth"],
       type = "l",
       col = "blue",
       xaxt = "n",
       xlab = "",
       xlim = c(10, 22),
       yaxt = "n",
       ylab = "",
       ylim = c(-24, 0),
       lwd = 4)
  # Add in axis for depth
  axis(2,
       at = seq(0, -24, by = -4),
       labels = seq(0, 24, by = 4))
  # Add in title for depth
  title(ylab = "Depth (m)",
        line = 1.5)
  # Add in axis for temperature
  axis(3,
       at = seq(10, 22, by = 2),
       labels = seq(10, 22, by = 2))
  mtext("Temperature (˚C)",
        line = 1.5,
        col = "blue")
  # Plot the DO values
  points(x = 10+WAVES.meta[which(WAVES.meta$Day == yearmonthday), "DO.(mg/L)"],
        y = -WAVES.meta[which(WAVES.meta$Day == yearmonthday), "Depth"],
        pch = 17)
  # Add axis for DO
  axis(1,
       at = seq(10, 22, by = 2),
       labels = seq(0, 12, by = 2))
  # Add title for DO
  title(x = "DO (mg/L)",
        line = 1.5)
  # Plot pH value
  points(x = WAVES.meta[which(WAVES.meta$Day == yearmonthday), "pH"]*15-108,
        y = -WAVES.meta[which(WAVES.meta$Day == yearmonthday), "Depth"],
        pch = 18,
        col = "red")
  # Axis for pH values
  axis(1,
       at = seq(10, 22, by = 2),
       labels = round((seq(10, 22, by = 2)+108)/15, 
                      2),
       line = 2.8,
       par(mgp = c(3, 0.5, 0)))
  #  Add title for pH 
  mtext("pH",
        line = 4.5,
        side = 1,
        col = "red")
  legend("bottomright",
         legend = c("Temperature",
                    "DO (mg/L)",
                    "pH"),
         col = c("blue",
                 "black",
                 "red"),
         pch = c(16, 17, 18),
         cex = 0.85)
  # Add title
#  title(main = paste("Lake Profile on ",
#                     "\n",
#                     substring(yearmonthday, 6, 7),
#                     "-",
#                     substring(yearmonthday, 9, 10),
#                     "-",
#                     substring(yearmonthday, 1, 4),
#                     sep = ""),
#        line = 3.3,
#        )
#  points(x = rep(16,
#                length(WAVES.meta[WAVES.meta$Day == yearmonthday, "Depth"])), 
#        y = -WAVES.meta[WAVES.meta$Day == yearmonthday, "Depth"],
#        pch = 16)
}

plotstuff("2009/09/02")
plotstuff("2009/09/30")
#####


# Generate a dendrogram of the samples
##########################
# Generate a distance matrix between the samples using the weighted data.
WAVES.REL.dist <- vegdist(t(WAVES.REL),
                          method = "bray")
WAVES.REL.dendro <- hclust(WAVES.REL.dist)
par(mar = c(2.5, 3, 1.8, 1))
plot(WAVES.REL.dendro,
     xlab = "",
     hang = -1,
     sub = "",
     main = "")
mtext("Sample Sites",
      side = 1,
      line = 1,
      cex = 1.2)
#####


# Generate a PCoA graph from the relativized data.
##########################
# Generate a distance matrix between the samples using the weighted data.
WAVES.REL.dist <- vegdist(t(WAVES.REL),
                          method = "bray")
# Generate the PCoA
WAVES.REL.pcoa <- cmdscale(WAVES.REL.dist)

# Set up the PCoA plot but, don't include the data yet.
plot(WAVES.REL.pcoa[, 1],
     WAVES.REL.pcoa[, 2],
     xlim = c(-0.25, 0.55),
     ylim = c(-0.35, 0.35),
     cex = 0,
     main = "WAVES data PCoA")

# Add data to the plot. 
text(WAVES.REL.pcoa[, 1],
     WAVES.REL.pcoa[, 2],
     c(2, 9, 10, 11, 12, 13, 16,
       2, 12, 14, 15, 16, 19),
     col = c(rep("blue", 7),
             rep("red", 6)),
     cex = 1)
#####


# Generate a phylum level stacked bar plot
##########################
# First things first. Make a new relative abundance table that uses
# the phylum-level classification rather than the OTU identity.
# Need to sum the rows with the same phylum classification.

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
melt.WAVES.REL.phyla[melt.WAVES.REL.phyla$value <= 0.005, "value"] <- 0
WAVES.REL.recast <- dcast(melt.WAVES.REL.phyla, 
                          formula = WAVES.phyla ~ variable)
WAVES.REL.recast <- WAVES.REL.recast[rowSums(WAVES.REL.recast[, -1]) != 0, ]
WAVES.REL.final <- melt(WAVES.REL.recast, 
                        id.vars=c("WAVES.phyla"))
# Remove unclassified group
WAVES.REL.final <- WAVES.REL.final[-grep("unclassified", WAVES.REL.final$WAVES.phyla), ]

# Phyla


mycolours <- c("#5276D8","#DF3B26","#F9A689","#7109AA","#61B7CF","#0C5AA6","#64DE89","#BDA71E","#E64319","#FF7A00","#5667AF","#403671")
phyla.figure <- ggplot(WAVES.REL.final, 
       aes(x=variable, 
           y=value, 
           fill=WAVES.phyla)) + 
  geom_bar(width=.5, colour = "black", stat="identity") +
  scale_y_continuous(labels = percent_format(),
                     breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                     name = "") + 
  scale_x_discrete(labels = c(2, 9, 10, 11, 12, 13, 16, 2, 12, 14, 15, 16, 19),
                   name = "Depth (m)") +
  theme(panel.background = element_rect(fill='white'), 
        axis.text.x = element_text(angle=0)) + 
  scale_fill_manual(values=mycolours)
ggsave(file = "plots/phyla_figure.eps")


# Split into A and B datas
# Day A
day.A <- WAVES.REL.final[1:77, ]
mycolours <- c("#5276D8","#DF3B26","#F9A689","#7109AA","#61B7CF","#0C5AA6","#64DE89","#BDA71E","#E64319","#FF7A00","#5667AF","#403671")
ggplot(day.A, 
       aes(x=variable, 
           y=value, 
           fill=WAVES.phyla)) + 
  geom_bar(width=.5, colour = "black", stat="identity") +
  scale_y_continuous(labels = percent_format(),
                     breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                     name = "") + 
  scale_x_discrete(labels = c(2, 9, 10, 11, 12, 13, 16),
                   name = "Depth (m)") +
  theme(panel.background = element_rect(fill='white'), 
        axis.text.x = element_text(angle=0)) + 
  scale_fill_manual(values=mycolours)
#+ theme(axis.title.x = element_text(vjust=2))

# Day B
day.B <- WAVES.REL.final[78:143, ]
mycolours <- c("#5276D8","#DF3B26","#F9A689","#7109AA","#61B7CF","#0C5AA6","#64DE89","#BDA71E","#E64319","#FF7A00","#5667AF","#403671")
ggplot(day.B, 
       aes(x=variable, 
           y=value, 
           fill=WAVES.phyla)) + 
  geom_bar(width=.5, colour = "black", stat="identity") +
  scale_y_continuous(labels = percent_format(),
                     breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                     name = "") + 
  scale_x_discrete(labels = c(2, 12, 14, 15, 16, 19),
                   name = "Depth (m)") +
  theme(panel.background = element_rect(fill='white'), 
        axis.text.x = element_text(angle=0)) + 
  scale_fill_manual(values=mycolours)
#+ theme(axis.title.x = element_text(vjust=2))


#####


# Generate a heat map with the top 10 orders from each sample
# Have a dendrogram at the top and bottom
##########################
# Generate an abundance table with all the phylogenetic classifications as the names
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
# Combine this file with the abundance table
WAVES.REL.tax <- cbind(WAVES.tax.list,
                       WAVES.REL)

# Write out a .csv file for this.
write.csv(WAVES.REL.tax,
          "dataEdited/WAVES.taxonony.abund",
          row.names = FALSE,
          quote = FALSE)

# Combine by order
WAVES.REL.species <- aggregate(cbind(A1002, A1009, A1010, A1011, A1012, A1013, A1216, B1002, B1012, B1014, B1015, B1016, B1019) ~ WAVES.REL.tax$Species,
                             data = WAVES.REL.tax,
                             FUN = sum)
# Rename the genus classification column and make into a character
names(WAVES.REL.species)[1] <- "Species"
WAVES.REL.species$Species <- as.character(WAVES.REL.species$Species)
rownames(WAVES.REL.species) <- WAVES.REL.species$Species
# Pull out the top 10 orders in each grouping.
top.ten <- matrix(0,
                  nrow = 10,
                  ncol = 13)
for (i in 2:ncol(WAVES.REL.species)) {
  top.ten[, (i-1)] <- WAVES.REL.species[head(sort(WAVES.REL.species[, i],
                                                decreasing = TRUE,
                                                index.return = TRUE)$ix,
                                           10), "Species"]
  }
# Name the columns for the samples
colnames(top.ten) <- site.names
# Find which ones are unique
top.ten <- unique(as.vector(top.ten))
# Remove unclassified samples
group.for.comparison <- top.ten[-grep("unclassified", top.ten)]
# Remove the "f__"
#for (i in 1:length(group.for.comparison)) {
#  group.for.comparison[i] <- 
#    print(strsplit(group.for.comparison[i],
#                                      "f__")[[1]])
#    }



# Make the heatmap
fig5B<-heatmap.2(apply(WAVES.REL.species[group.for.comparison, -1],
                       1,
                       sqrt),
                 col=hc(20),
                 #scale="col",
                 key=TRUE,
                 symkey=FALSE,
                 trace="none",
                 density.info="none",
                 dendrogram="both",
                 labCol = group.for.comparison,
                 #colCol=c(map$Class_Color),
                 margins=c(5,13), srtCol=90)



# Let's try this with a shared table.

WAVES.shared <- read.table("dataEdited/WAVES.final.shared",
                           # OTU numbers are the heading names
                           header = TRUE,
                           # The file is tab-delimited
                           sep = "\t",
                           # Use the site names as row names
                           row.names = 2)
WAVES.shared <- WAVES.shared[, -c(1,2)]











# The following was commented out because it is not an appropriate test.
# All of the factors we're looking at are correlated to depth. Since the 
# variables are not independent of that gradient, this is not appropropriate.

# Correlate the environmental data to the PCoA components
#WAVES.env <- cbind(round(cor(WAVES.meta[, -1],
#                             WAVES.REL.pcoa[, 1]),
#                         2),
#                   round(cor(WAVES.meta[, -1],
#                             WAVES.REL.pcoa[, 2]),
#                         2))

# Fit the environmental factors onto the ordination. 
#WAVESEF <- envfit(WAVES.REL.pcoa,
#                  WAVES.meta[, -1])
# Add that fit to the plot
#plot(WAVESEF)
#####





# Check out phylum-level variation in community composition
##########################
# First, clean up taxonomy file
#test <- WAVES.tax.phylo[c(1:10), ]
#WAVES.tax.phylo <- 
# Need to make a phyloseq object.
# First, turn them into a tax_table and otu_table
#WAVES.tax.phylo <- tax_table(WAVES.tax)
#WAVES.REL.phylo <- otu_table(WAVES.REL, taxa_are_rows = TRUE)
#WAVES.phyloseq <- phyloseq(WAVES.REL.phylo,
#                           WAVES.tax.phylo)
#plot_bar(WAVES.phyloseq)
#
#head(WAVES.tax.phylo)

# Shared table
#WAVES.shared <- read.table("dataEdited/WAVES.final.shared",
#                           # OTU numbers are the heading names
#                           header = TRUE,
#                           # The file is tab-delimited
#                           sep = "\t",
#                           # Use the site names as row names
#                           row.names = 2)
#WAVES.shared <- WAVES.shared[, -c(1,2)]
# Transpose matrix
#WAVES.shared <- t(WAVES.shared)
