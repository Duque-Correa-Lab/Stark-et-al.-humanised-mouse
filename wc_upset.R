# Load required libraries
library(ComplexUpset)
library(reshape2)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(ggprism)
library(writexl)


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

# Convert abundance to presence/absence (1/0)
t.species.df = as.data.frame(t(species.df))
t.species.df[t.species.df > 0] = 1 # filter for relative abundance

# Add sample annotation and aggregate
t.species.df$Annot = metadata.filt[rownames(t.species.df), "sample.title"]
annot.agg = aggregate(. ~ Annot, data=t.species.df, FUN=sum)
rownames(annot.agg) = annot.agg[,1]
annot.agg = annot.agg[,-1]

# Transpose and filter for non-zero rows
upset.df = t(annot.agg)
upset.df[upset.df > 0] = 1
upset.df = as.data.frame(upset.df[which(rowSums(upset.df > 0) > 0),])

# Plot basic upset plot (set size + intersection matrix)
upset.plot = upset(upset.df, colnames(upset.df), name="")

# Display plot
upset.plot

# Add taxonomic ranks
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
upset.df$Lineage = rownames(upset.df)
upset.tax = separate(data = upset.df, col = Lineage, sep = ";", into = ranks)

# Get top 12 families by frequency
upset.tax.sums <- table(upset.tax$Family)
upset.tax.sums.df <- as.data.frame(upset.tax.sums)
top.family = upset.tax.sums.df[order(upset.tax.sums.df$Freq, decreasing=TRUE)[1:12],"Var1"]
top.family_df = rownames(upset.tax[upset.tax$Family %in% top.family,])
top.family_df <- upset.tax[upset.tax$Family %in% top.family, ]

# Clean up taxonomic prefixes
for (tax in ranks){
  top.family_df[,tax] = gsub("[dpcofgs]__", "", top.family_df[,tax])
}

top.family_df = top.family_df[,c("WT", "D2 HMA","D7 HMA","Domain", "Phylum","Class","Order","Family","Genus","Species")]

# plot upset tax
paired_colors <- brewer.pal(12, "Paired")
names(paired_colors) <- unique(top.family_df$Family)
paired_colors[1] <- "#E58606FF"
paired_colors[2] <- "#2C3E50"
paired_colors[3] <- "#52BCA3FF"
paired_colors[4] <- "#99C945FF"
paired_colors[5] <- "#CC61B0FF"
paired_colors[6] <- "#24796CFF"
paired_colors[7] <- "#DAA51BFF"
paired_colors[8] <- "#2F8AC4FF"
paired_colors[9] <- "#764E9FFF"
paired_colors[10] <- "#ED645AFF"
paired_colors[11] <- "#CC3A8EFF"
paired_colors[12] <- "#A5AA99FF"


# --------------------------------------------------
# Create final UpSet plot with taxonomic annotations
# --------------------------------------------------

# Create barplot annotation using ComplexUpset-compatible format
upset.plot.tax.top = upset(top.family_df, c("WT", "D2 HMA", "D7 HMA"), name="",
                           base_annotations=list(), 
                           annotations = list(
                             taxon=(
                               ggplot(mapping=aes(fill=Family))
                               + geom_bar(stat='count', position='fill', colour="black", linewidth=0.1)
                               + scale_y_continuous(labels=scales::percent_format())
                               + scale_fill_manual(values = paired_colors, name="Family")
                               + theme(axis.title.y = element_text(size = 16))
                               + theme(axis.text.y = element_text(size = 13), face = "bold")
                               + labs(x = "Number of species", y = "")
                               + theme_prism()
                               + theme(legend.title = element_text(size = 20, face = "bold"))
                               + theme (axis.title.x = element_blank(),
                                        legend.position = "left", legend.text = element_text(size = 18, face = "italic")))),
                           width_ratio = 0.25,
                           themes=upset_default_themes(text=element_text(size=12)),
                           stripes = upset_stripes(colors = c("indianred3","steelblue3", "gray43"))
)

# Display plot
upset.plot.tax.top

