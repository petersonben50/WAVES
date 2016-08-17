# Always get a clean start
rm(list = ls())

# Change working directory
##########################
setwd("/Users/benjaminpeterson/Documents/gradSchool/research/mercury/WAVES")
#####

# Load up needed packages
##########################
library("vegan")
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

# Read in OTU classification table
WAVES.otus <- read.csv2("dataEdited/WAVES.taxonomy",
                        stringsAsFactors = FALSE,
                        header = FALSE)
# Make OTU names the row names
otu.names <- WAVES.otus[, 1]
rownames(WAVES.otus) <- otu.names
WAVES.otus <- WAVES.otus[, -1]

# Read in the .csv file with the environmental data
# Transpose it so that the samples are the rows.
WAVES.meta <- as.data.frame(t(read.csv("dataRaw/2,30 Sept 2009 ME canoco.csv",
                       header = TRUE,
                       row.names = 1,
                       stringsAsFactors = FALSE)))
# Keep only the values for the samples that we sequenced.
WAVES.meta <- WAVES.meta[row.names(WAVES.REL), ]
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

# Overview of Environmental Data
##########################
plotstuff <- function(yearmonthday) {
  # Set up plot parameters
  par(mar = c(6, 3, 5, 2),
      mgp = c(3, 0.5, 0))
  
  # Set up plot, put temperature data points on it.
  plot(WAVES.meta[which(WAVES.meta$Day == yearmonthday), "Temp"],
       -WAVES.meta[which(WAVES.meta$Day == yearmonthday), "Depth"],
       type = "l",
       col = "blue",
       xaxt = "n",
       xlab = "",
       xlim = c(10, 22),
       yaxt = "n",
       ylab = "",
       ylim = c(-20, 0),
       lwd = 4)
  # Add in axis for depth
  axis(2,
       at = seq(0, -20, by = -4),
       labels = seq(0, 20, by = 4))
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
  lines(x = 10+WAVES.meta[which(WAVES.meta$Day == yearmonthday), "DO.(mg/L)"],
        y = -WAVES.meta[which(WAVES.meta$Day == yearmonthday), "Depth"],
        lwd = 4)
  # Add axis for DO
  axis(1,
       at = seq(10, 22, by = 2),
       labels = seq(0, 12, by = 2))
  # Add title for DO
  title(x = "DO (mg/L)",
        line = 1.5)
  # Plot pH value
  lines(x = WAVES.meta[which(WAVES.meta$Day == yearmonthday), "pH"]*15-108,
        y = -WAVES.meta[which(WAVES.meta$Day == yearmonthday), "Depth"],
        lwd = 4,
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
  # Add title
  title(main = paste("Lake Profile on ",
                     substring(yearmonthday, 6, 7),
                     "-",
                     substring(yearmonthday, 9, 10),
                     "-",
                     substring(yearmonthday, 1, 4),
                     sep = ""),
        line = 3.3)
}

plotstuff("2009/09/02")
plotstuff("2009/09/30")
#####

# Generate a PCoA graph from the relativized data.
##########################
# Generate a distance matrix between the samples using the weighted data.
WAVES.REL.dist <- vegdist(WAVES.REL,
                          method = "bray")
# Generate the PCoA
WAVES.REL.pcoa <- cmdscale(WAVES.REL.dist)

# Set up the PCoA plot but, don't include the data yet.
plot(WAVES.REL.pcoa[, 1],
     WAVES.REL.pcoa[, 2],
     cex = 0,
     main = "WAVES data PCoA")

# Add data to the plot. 
text(WAVES.REL.pcoa[, 1],
     WAVES.REL.pcoa[, 2],
     attr(WAVES.REL.pcoa, "dimnames")[[1]],
     cex = 0.6)

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
