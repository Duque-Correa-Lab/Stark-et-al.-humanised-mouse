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
setwd("/your/path/here/)
metadata = read.delim("MSvsLC_metadata.txt")
rownames(metadata) = metadata$run.accession

# Load abundance data
in_file = "merged_abundance_table.txt"
data_otu = read.delim(in_file, skip=1L)

# Keep only species-level entries (exclude taxonomy summary rows)
species.df = data_otu[grepl("s__", data_otu$clade_name),]
species.df = species.df[!grepl("t__", species.df$clade_name),]
rownames(species.df) = species.df$clade_name
species.df = species.df[,-1]

# Subset to D7 HMA samples only
common.samples = intersect(colnames(species.df), rownames(metadata))
metadata.subset = metadata[which(metadata$annot1 != "WT" & metadata$run.accession %in% common.samples),]
metadata.subset.D7HMA <- metadata.subset[metadata.subset$annot1 != "D2 HMA", ]
species.subset = species.df[,rownames(metadata.subset)]

# ----------------------------------------
# Run Maaslin2
# ----------------------------------------

# Output directory for Maaslin2 results
output_dir = "maaslin2_MSvsLC"

# Run MaAsLin2 model with specified options
cat("MaAsLin2...\n")
fit_data = Maaslin2(cores = 4, normalization = "NONE", transform = "LOG",
                    min_prevalence = 0.01, max_significance = 0.05, correction = "BH",
                    input_data = species.subset, input_metadata = metadata.subset, output = output_dir,
                    fixed_effects = c("annot2"),
                    reference = c("annot2,mucosal scraping"),
                    plot_heatmap = TRUE, plot_scatter = TRUE)

# ----------------------------------------
# Analysis and Plotting for target species
# ----------------------------------------

results_dir <- "/Users/starklaraliz/Library/CloudStorage/OneDrive-UniversityofCambridge/untitled_folder/GB/Cambridge/PhD/PhD2/HMA project/Metagenomics/clean code/MSvsLC/maaslin2_MSvsLC"
norm_data = t(read.table(paste0(results_dir, "/features/filtered_data_norm.tsv")))
colnames(norm_data) = norm_data[1,]
rownames(norm_data) = norm_data[,1]
tax_list = rownames(norm_data)[2:length(rownames(norm_data))]
norm_data = apply(norm_data[2:length(rownames(norm_data)),2:length(colnames(norm_data))], 2, as.numeric)
norm_data = as.data.frame(norm_data)
rownames(norm_data) = tax_list
norm_data$Species = sub(".+\\.s__(.*)", "\\1", rownames(norm_data))
norm_data$Species = gsub("\\.", " ", norm_data$Species)
target_sp = c("Acetatifactor sp910579755", "COE1 sp910579275", "UBA3282 sp910584965", "Anaerotruncus sp000403395")


for (sp in target_sp){
  plotdf = melt(norm_data[norm_data$Species == sp,])
  plotdf$annot2 = metadata[plotdf$variable, "annot2"]
  plotdf$annot2 = factor(plotdf$annot2, levels = c("luminal content", "mucosal scraping"))
   # Create summary for plot
  summary_df <- plotdf %>%
    group_by(annot2) %>%
    summarise(
      mean_value = mean(plotdf$value),
      standard_error = sd(plotdf$value)/sqrt(length(plotdf$value)),
      ymin = mean_value - standard_error,
      ymax = mean_value + standard_error,
      Q1 = mean_value - standard_error,
      Q3 = mean_value + standard_error
    )

  p = ggplot(plotdf, aes(annot2, value, color = annot2)) +
    
    # Add mean + SE
    stat_summary(fun = mean, geom = "crossbar", color = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, size = 1, color = "black") +
    
    # Custom theme and label formatting
    theme_prism() +
    labs(title = "Species Comparison", x = "", y = "Relative abundance (%)") +
    
    # Color dots by group
    geom_jitter(height = 0, width = 0.2, size = 2) +  # Add size for clarity
    
    scale_color_manual(values = c("luminal content" = "#FA8", "mucosal scraping" = "#9E2A2F")) +
    
    # Axis and text styling
    theme(axis.title = element_text(size=14, face = "bold")) + 
    theme(axis.text = element_text(size=14, face = "bold")) +
    
    theme(legend.position = "none") +
    theme(plot.title = element_blank()) +
    
    scale_y_continuous(n.breaks = 8, labels = scales::number_format(accuracy = 0.01), 
                       expand = expansion(mult = c(0.05, 0.1))) +
    
    stat_compare_means(comparisons = list(c("luminal content", "mucosal scraping")),
                       method = "t.test",
                       tip.length = 0,
                       map_signif_level = FALSE,
                       bracket.size = 1,
                       size = 6)
  
  print(p)
  }



