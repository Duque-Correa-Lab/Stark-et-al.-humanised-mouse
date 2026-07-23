# Load required libraries
library(Maaslin2)   # For differential abundance testing
library(reshape2)   # For data reshaping (melt function)
library(ggplot2)    # For visualization
library(ggprism)    # For prism-style ggplots
library(dplyr)      # For data manipulation
library(ggpubr)     # For stat_compare_means()

# ----------------------------------------
# Load and preprocess data
# ----------------------------------------

# Set working directory and load metadata
cat("Loading data ...\n")
setwd("/your/path/here/")
metadata = read.delim("metadata_WC_HMAvsWT.txt")
rownames(metadata) = metadata$run.accession

# Load abundance data
in_file = "merged_abundance_table.txt"
data_otu = read.delim(in_file, skip=1L)

# Keep only species-level entries (exclude taxonomy summary rows)
species.df = data_otu[grepl("s__", data_otu$clade_name),]
species.df = species.df[!grepl("t__", species.df$clade_name),]
rownames(species.df) = species.df$clade_name
species.df = species.df[,-1]

# ----------------------------------------
# Run Maaslin2
# ----------------------------------------

# Output directory for Maaslin2 results
output_dir = "maaslin2_WC_HMAvsWT"

# Run MaAsLin2 model with specified options
fit_data = Maaslin2(cores = 4, normalization = "NONE", transform = "LOG",
                    min_prevalence = 0.01, max_significance = 0.05, correction = "BH",
                    input_data = species.df, input_metadata = metadata, output = output_dir,
                    fixed_effects = c("sample.title"),
                    reference = c("sample.title,D7 HMA"),
                    plot_heatmap = TRUE, plot_scatter = TRUE)

# ----------------------------------------
# Analysis and Plotting for target species
# ----------------------------------------

results_dir <- "/Users/starklaraliz/Library/CloudStorage/OneDrive-UniversityofCambridge/untitled_folder/GB/Cambridge/PhD/PhD2/HMA project/Metagenomics/clean code/HMAvsWT/maaslin2_WC_HMAvsWT"
norm_data = t(read.table(paste0(results_dir, "/features/filtered_data_norm.tsv")))
colnames(norm_data) = norm_data[1,]
rownames(norm_data) = norm_data[,1]
tax_list = rownames(norm_data)[2:length(rownames(norm_data))]
norm_data = apply(norm_data[2:length(rownames(norm_data)),2:length(colnames(norm_data))], 2, as.numeric)
norm_data = as.data.frame(norm_data)
rownames(norm_data) = tax_list
norm_data$Species = sub(".+\\.s__(.*)", "\\1", rownames(norm_data))
norm_data$Species = gsub("\\.", " ", norm_data$Species)
target_sp = c("Bacteroides clarus", "Bacteroides xylanisolvens", "MGBC163490 sp910588075", "COE1 sp910579275", "Bacteroides cellulosilyticus")

for (sp in target_sp){
  plotdf = melt(norm_data[norm_data$Species == sp,])
  plotdf$sample.title = metadata[plotdf$variable, "sample.title"]
  plotdf$sample.title = factor(plotdf$sample.title, levels = c("WT", "D2 HMA", "D7 HMA"))
  summary_df <- plotdf %>%
    group_by(sample.title) %>%
    summarise(
      mean_value = mean(plotdf$value),
      standard_error = sd(plotdf$value)/sqrt(length(plotdf$value)),
      ymin = mean_value - standard_error,
      ymax = mean_value + standard_error,
      Q1 = mean_value - standard_error,
      Q3 = mean_value + standard_error
    )
  
  p <- ggplot(plotdf, aes(x = sample.title, y = value, fill = sample.title)) +
    
    # Remove default boxplot and manually add mean + SE
    stat_summary(fun = mean, geom = "crossbar") +  # Mean point
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, size = 1, color = "black") +  # SE error bars
    
    # Optional jitter plot to show raw data
    geom_jitter(aes(colour = sample.title), height = 0, width = 0.2, alpha = 1, size = 2) +
    
    # Custom styling
    theme_prism() +
    labs(title = unique(plotdf$Species), x = "", y = "Relative abundance (%)") +
    scale_color_manual("sample.title", values = c("WT" = "gray43", "D2 HMA" = "steelblue3", "D7 HMA" = "indianred3")) + 

    # Theme modifications
    theme(axis.title = element_text(size = 14, face = "bold")) +
    theme(axis.text = element_text(size = 14, face = "bold")) +
    theme(plot.title = element_blank()) +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01))
  
  print(p)
  }

