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
setwd("/your/path/here/")

in_file = "merged_abundance_table.txt"
data_otu = read.delim(in_file, skip=1L)

# Filter OTU table to species-level entries
species.df = data_otu[grepl("s__", data_otu$clade_name),]
species.df = species.df[!grepl("t__", species.df$clade_name),]
rownames(species.df) = species.df$clade_name
species.df = species.df[,-1]

# Load metadata and align with filtered OTU data
metadata = read.delim("metadata_WC_HMAvsWT.txt")
rownames(metadata) = metadata$run.accession
metadata.filt = metadata[colnames(species.df),]

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
betadiv.pca = data.frame(row.names=rownames(pcoa.axes), # add metadata columns here
                         PC1=pcoa.axes[,1],
                         PC2=pcoa.axes[,2],
                         Annot=metadata.filt[rownames(pcoa.axes),"sample.title"])

# ----------------------------------------
# Plotting PCoA
# ----------------------------------------

pca.plot = ggplot(betadiv.pca, aes(x=PC1, y=PC2, colour=Annot)) + 
  geom_point(size=3, alpha=1) +
  theme_classic() +
  theme_prism() +
  theme(panel.grid.minor = element_blank()) + 
  ylab(paste("PC2"," (",pc2_var,"%)",sep="")) + 
  xlab(paste("PC1"," (",pc1_var,"%)",sep="")) +
  theme(legend.position="top", legend.text=element_text(size=18, face = "bold"),
        legend.title = element_blank()) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(axis.title = element_text(size=18, face = "bold")) + 
  theme(axis.text = element_text(size=14, face = "bold")) + 
  scale_colour_manual(values = c("WT" = "gray43", "D2 HMA" = "steelblue3", "D7 HMA" = "indianred3"),
                      labels = c("WT", "D2 HMA", "D7 HMA"),
                      breaks = c("WT", "D2 HMA", "D7 HMA"))

# Display plot
pca.plot


