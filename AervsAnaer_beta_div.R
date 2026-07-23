# Load required libraries
library(ape)            # For PCoA (principal coordinates analysis)
library(vegan)          # For ecological analysis (e.g., Bray-Curtis distances)
library(ggplot2)        # Plotting
library(tidyr)          # Tidy data handling
library(data.table)     # Efficient data handling
library(ggprism)        # Publication-ready ggplot themes

# ----------------------------------------
# Load and preprocess data
# ----------------------------------------

# Set working directory and load OTU data
cat("Loading data ...\n")
setwd("/your/path/here/") # <-- Use your actual path
data_otu = read.delim("merged_abundance_table.txt", skip=1L)

# Filter OTU table to species-level entries
species.df = data_otu[grepl("s__", data_otu$clade_name),]
species.df = species.df[!grepl("t__", species.df$clade_name),]
rownames(species.df) = species.df$clade_name
species.df = species.df[,-1]

# Load metadata and align with filtered OTU data
metadata = read.delim("AerobevsAnaerobe_scraping_metadata.txt")
rownames(metadata) = metadata$run.accession
metadata.filt = metadata[colnames(species.df),]
metadata.filt = metadata.filt[metadata.filt$annot1 != "WT", ]
metadata.filt = metadata.filt[metadata.filt$annot1 != "D2 HMA", ]

# ----------------------------------------
# Beta diversity analysis (Bray-Curtis + PCoA)
# ----------------------------------------

# Calculate Bray-Curtis distance between samples
cat("Calculating distance matrices ...\n")
dist.df = vegdist(t(species.df), method="bray")

# Perform Principal Coordinates Analysis (PCoA)
cat("Building PCoA ...\n")
pcoa.data = pcoa(dist.df)
pcoa.axes = data.frame(pcoa.data$vectors)
pc1_var = round(pcoa.data$values[1,2]*100,1)
pc2_var = round(pcoa.data$values[2,2]*100,1)

# Merge metadata with PCoA results
betadiv.pca = data.frame(row.names=rownames(pcoa.axes), # add metadata columns here
                         PC1=pcoa.axes[,1],
                         PC2=pcoa.axes[,2],
                         Annot1=metadata.filt[rownames(pcoa.axes),"annot1"],
                         Annot2=metadata.filt[rownames(pcoa.axes),"annot2"])

# Remove samples with missing annotation
betadiv.pca_clean <- betadiv.pca[!is.na(betadiv.pca$Annot2), ]

# ----------------------------------------
# Plotting PCoA
# ----------------------------------------

pca.plot = ggplot(betadiv.pca_clean, aes(x=PC1, y=PC2, colour=Annot1, shape=Annot2)) + 
  geom_point(size=3, alpha=0.4) +
  labs(shape = "AerobicvsAnaerobic") +
  theme_prism() +
  theme(panel.grid.minor = element_blank()) + 
  ylab(paste("PC2"," (",pc2_var,"%)",sep="")) + 
  xlab(paste("PC1"," (",pc1_var,"%)",sep="")) +
  theme(legend.position="right", legend.text=element_text(size=14, face = "bold"),
        legend.title = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.title.y = element_text(size=18)) + 
  theme(axis.text.y = element_text(size=14)) + 
  theme(axis.title.x = element_text(size=18)) + 
  theme(axis.text.x = element_text(size=14)) +
  scale_colour_manual(values = c("D7 HMA" = "indianred4")) +
  stat_ellipse(size = 0.5, level = 0.95, linetype = "dashed")

# Display plot
pca.plot

