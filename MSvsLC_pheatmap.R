# Load required libraries
library(pheatmap)     # For heatmap plotting
library(stringr)      # For string manipulation
library(grid)         # For plot layout (not directly used but required by pheatmap)
library(ggplot2)      # For color tools (e.g., alpha transparency)

# ----------------------------------------
# Load and preprocess data
# ----------------------------------------

# Set working directory and load metadata
cat("Loading data ...\n")
setwd("/your/path/here/")
significant_otu = read.delim("maaslin2_MSvsLC/significant_results.tsv")
significant_otu$feature = gsub("\\.([pcofgs]__)", ";\\1", significant_otu$feature)
significant_otu$feature = gsub("\\.", " ", significant_otu$feature)

# Load full species abundance data and metadata
data_otu = read.delim("merged_abundance_table.txt", skip=1L)

metadata = read.delim("MSvsLC_metadata.txt")
metadata = metadata[!(metadata$annot1 %in% c("D2 HMA", "WT")), ]


# ----------------------------------------
# Prepare species abundance data for heatmap
# ----------------------------------------

# Keep only species that were significant in Maaslin2
data_otu <- data_otu[data_otu$clade_name %in% significant_otu$feature, ]

#draw the heatmap
rownames(data_otu) = data_otu$clade_name
data_otu = data_otu[,-1]
data_otu <- data_otu[,metadata$run.accession]

# Sample-level annotation (column annotation in heatmap)
sample_coloring = data.frame(Group = as.vector(metadata$annot1),
                             sample_name = metadata$run.accession,
                             Source = metadata$annot2)

rownames(sample_coloring) <- metadata$run.accession
sample_coloring <- sample_coloring[,"Source", drop = FALSE]

# Keep only top 10 most significant species
data_otu_first10 <- head(data_otu, 10)

# Extract taxonomy from rownames
specm = str_split_fixed(rownames(data_otu_first10), ";", 7)
specdf = as.data.frame(gsub("[dpcofgs]__", "", specm))
colnames(specdf) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(specdf) <- rownames(data_otu_first10)

family_df <- as.data.frame(specdf[,"Family"])
colnames(family_df) <- c("Family")
rownames(family_df) <- rownames(data_otu_first10)

# Semi-transparent colors
transparent_content <- alpha("#FA8", 0.8)
transparent_scraping <- alpha("#9E2A2F", 0.8)

# Okabe-Ito colorblind-safe palette
palette_OkabeIto <- c("#E58606FF", "#2C3E50", "#52BCA3FF", "#99C945FF",
                      "#CC61B0FF",  "#24796CFF", "#DAA51BFF", "#2F8AC4FF", "#764E9FFF",
                      "#ED645AFF", "#CC3A8EFF", "#A5AA99FF")
palette_OkabeIIto <- c(transparent_content, transparent_scraping)

# Assign colors directly without interpolation
ann_colors <- list(
  Source = setNames(palette_OkabeIIto[seq_along(unique(sample_coloring$Source))], unique(sample_coloring$Source)),
  Family = setNames(palette_OkabeIto[seq_along(unique(family_df$Family))], unique(family_df$Family))
)

names(ann_colors$Source) <- unique(sample_coloring$Source)
names(ann_colors$Family) <- unique(family_df$Family)

# ----------------------------------------
# Plot heatmap
# ----------------------------------------

p <- pheatmap(data_otu_first10,
                    scale = "row",
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    annotation_col = sample_coloring,
                    annotation_names_col = FALSE,
                    annotation_row = family_df,
                    annotation_names_row = FALSE,
                    annotation_colors = ann_colors,
                    labels_row = specdf$Species,
                    fontsize_row = 10,
                    cellwidth = 10,  # Set cell width
                    cellheight = 15, # Set cell height
)

# Display plot
p




