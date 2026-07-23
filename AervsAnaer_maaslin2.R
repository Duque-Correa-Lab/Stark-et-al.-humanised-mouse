# Load required libraries
library(Maaslin2)    # For multivariable association discovery
library(metagMisc)   # Misc tools for microbiome analysis

# ----------------------------------------
# Load and preprocess data
# ----------------------------------------

# Set working directory and load metadata
cat("Loading data ...\n")
setwd("/your/path/here/") # <-- Use your actual path
metadata = read.delim("AerobevsAnaerobe_scraping_metadata.txt")
rownames(metadata) = metadata$run.accession

# Load and filter species-level abundance data
in_file = "merged_abundance_table.txt"
data_otu = read.delim(in_file, skip=1L)
# colnames(data_otu) = gsub("_results_gtdb", "", colnames(data_otu))

# Keep only species-level entries (exclude taxonomy summary rows)
species.df = data_otu[grepl("s__", data_otu$clade_name),]
species.df = species.df[!grepl("t__", species.df$clade_name),]
rownames(species.df) = species.df$clade_name
species.df = species.df[,-1]

# Subset to D7 HMA samples only
common.samples = intersect(colnames(species.df), rownames(metadata))
metadata.subset = metadata[which(metadata$annot1 != "WT" & metadata$run.accession %in% common.samples),]
metadata.subset.D7HMA <- metadata.subset[metadata.subset$annot1 != "D2 HMA", ]
species.subset = species.df[,rownames(metadata.subset.D7HMA)]

# ----------------------------------------
# Run Maaslin2
# ----------------------------------------

# Output directory for Maaslin2 results
output_dir = "maaslin2_AervsAnaer_D7HMA_uncultured"

# Run MaAsLin2 model with specified options
cat("MaAsLin2...\n")
fit_data <- Maaslin2(
  input_data = species.subset,
  input_metadata = metadata.subset.D7HMA,
  output = output_dir,
  fixed_effects = c("annot2"),               # Test aerobic vs anaerobic (categorical)
  reference = c("annot2,Uncultured"),        # Set Uncultured as the reference level
  normalization = "NONE",                    # Assume input data already normalized
  transform = "LOG",                         # Apply log transformation
  min_prevalence = 0.01,                     # Filter features by prevalence
  max_significance = 0.05,                   # Adjusted p-value cutoff for results
  correction = "BH",                         # Use Benjamini-Hochberg correction
  plot_heatmap = TRUE,                       # Generate heatmap for significant features
  plot_scatter = TRUE,                       # Generate scatter plots for top associations
  cores = 4                                  # Use 4 CPU cores
)





