# Load required libraries
library(pheatmap)     # For heatmap plotting
library(stringr)      # For string manipulation
library(grid)         # For plot layout (not directly used but required by pheatmap)
library(ggplot2)      # For color tools (e.g., alpha transparency)

# ----------------------------------------
# Load and preprocess data
# ----------------------------------------

# Set working directory to where input files are located
setwd("/your/path/here/") # <-- Use your actual path

# Load significant species identified by Maaslin2
cat("Loading data ...\n")
significant_otu = read.delim("maaslin2_AervsAnaer_D7HMA_anaerobe/significant_results.tsv")
significant_otu$feature = gsub("\\.([pcofgs]__)", ";\\1", significant_otu$feature)
significant_otu$feature = gsub("\\.", " ", significant_otu$feature)

# Load full species abundance data and metadata
data_otu = read.delim("merged_abundance_table.txt", skip=1L)

metadata = read.delim("AerobevsAnaerobe_scraping_metadata.txt")
metadata = metadata[!(metadata$annot1 %in% c("D2 HMA", "WT")), ]

# ----------------------------------------
# Prepare species abundance data for heatmap
# ----------------------------------------

# Keep only species that were significant in Maaslin2
data_otu <- data_otu[data_otu$clade_name %in% significant_otu$feature, ]

# Set species names as rownames
rownames(data_otu) = data_otu$clade_name
data_otu = data_otu[,-1]

# Subset to match only selected samples in metadata
data_otu <- data_otu[,metadata$run.accession]

# Sample-level annotation (column annotation in heatmap)

sample_coloring = data.frame(Group = as.vector(metadata$annot1),
                             sample_name = metadata$run.accession,
                             Source = metadata$annot2)
rownames(sample_coloring) <- metadata$run.accession
sample_coloring <- sample_coloring[,"Source", drop = FALSE]

# Keep only top 30 most significant species
data_otu_first30 <- head(data_otu, 30)

# Extract taxonomy from rownames
specm = str_split_fixed(rownames(data_otu_first30), ";", 7)
specdf = as.data.frame(gsub("[dpcofgs]__", "", specm))
colnames(specdf) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(specdf) <- rownames(data_otu_first30)

family_df <- as.data.frame(specdf[,"Family"])
colnames(family_df) <- c("Family")
rownames(family_df) <- rownames(data_otu_first30)

# Semi-transparent colors
transparent_gray12 <- alpha("gray12", 0.5)
transparent_indianred4 <- alpha("#8B3A3A", 0.5)

# Okabe-Ito colorblind-safe palette
palette_OkabeIto <- c("#E58606FF", "#2C3E50", "#52BCA3FF", "#99C945FF",
                      "#CC61B0FF",  "#24796CFF", "#DAA51BFF", "#2F8AC4FF", "#764E9FFF",
                      "#ED645AFF", "#CC3A8EFF", "#A5AA99FF")
palette_OkabeIIto <- c(transparent_indianred4, "white", transparent_gray12)

# Assign colors directly without interpolation
ann_colors <- list(
  Source = setNames(palette_OkabeIIto[seq_along(unique(sample_coloring$Source))], unique(sample_coloring$Source)),
  Family = setNames(palette_OkabeIto[seq_along(unique(family_df$Family))], unique(family_df$Family))
)

names(ann_colors$Source) <- unique(sample_coloring$Source)
names(ann_colors$Family) <- unique(family_df$Family)

# Custom sample order
order <- c("ERR12418660", "ERR12418661", "ERR12418667", "ERR12418655", "ERR12418670", "ERR12418672", "ERR12418663", "ERR12418665", "ERR12418658", "ERR12418659", 
           "ERR12418601","ERR12418605","ERR12418622", "ERR12418626", "ERR12418629", "ERR12418631", "ERR12418608", "ERR12418610", "ERR12418592", "ERR12418594",
           "ERR12418602", "ERR12418606", "ERR12418623", "ERR12418627", "ERR12418633", "ERR12418635", "ERR12418612", "ERR12418614", "ERR12418596", "ERR12418598")
df = data_otu_first30[order]

# ----------------------------------------
# Plot heatmap
# ----------------------------------------

p <- pheatmap(
  log10(df + 1),                     # Log-transform data for visualization
  show_rownames = TRUE,              # Show species names
  show_colnames = FALSE,             # Hide sample names
  cluster_cols = TRUE,               # Enable clustering of columns
  color = colorRampPalette(c("white", "red"))(100),  # Heatmap color gradient
  annotation_col = sample_coloring,  # Sample group annotation
  annotation_names_col = FALSE,
  annotation_row = family_df,        # Family-level row annotation
  annotation_colors = ann_colors,    # Custom annotation colors
  annotation_names_row = FALSE,
  labels_row = specdf$Species,       # Label rows with species names
  fontsize_row = 10,
  cellwidth = 10,
  cellheight = 15
)

# Display plot
p
