# Load required libraries
library(vegan)         # For diversity indices (e.g., Shannon)
library(ggplot2)       # For plotting
library(stringr)       # String manipulation
library(dunn.test)     # For post-hoc pairwise testing
library(dplyr)         # Data manipulation
library(ggprism)       # For prism-style ggplots
library(ggpubr)        # For geom_signif()

# ----------------------------------------
# Load and preprocess species-level OTU data
# ----------------------------------------

# Set working directory and load OTU data
cat("Loading data ...\n")
setwd("your path")
data_otu = read.delim("merged_abundance_table_whole_caecum.txt", skip=1L)

# Extract and clean taxonomy information
phylum.df = as.data.frame(data_otu[grepl("s__", data_otu$clade_name), "clade_name"])
colnames(phylum.df) = c("clade_name")
phylum.df = phylum.df[!grepl("t__", phylum.df$clade_name),]
phylum.df = str_split_fixed(phylum.df, ";", 7)
phylum.df = as.data.frame(sub("[dpcofgs]__", "", phylum.df))
colnames(phylum.df) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
phylum.df[length(rownames(phylum.df))+1,] = rep("Other", 7)

# Filter OTU table to species-level entries
species.df = data_otu[grepl("s__", data_otu$clade_name),]
species.df = species.df[!grepl("t__", species.df$clade_name),]
rownames(species.df) = species.df$clade_name
species.df = t(species.df[,-1])
rownames(species.df) <- sub("_results_gtdb$", "", rownames(species.df))

# Load and align metadata with species data
metadata = read.delim("metadata_WC_HMAvsWT.txt")
rownames(metadata) = metadata$run_accession
metadata.filt = metadata[rownames(species.df), ]


# ----------------------------------------
# Alpha diversity analysis (Shannon-index)
# ----------------------------------------

# Calculate Shannon diversity and merge with metadata
cat("Calculating Shannon-index ...\n")
alpha.div = diversity(species.df, index = "shannon")
metadata.filt$Shannon = alpha.div[rownames(metadata.filt)]
# metadata.filt$annot = factor(metadata.filt$sample_title, levels = c("WT", "D2 HMA", "D7 HMA"))

# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Shannon ~ annot, data = metadata.filt)

# Post-hoc Dunn’s test and adjust p-values if Kruskal-Wallis is significant
if (kruskal_test_result$p.value < 0.05) {
  dunn_test <- dunn.test(metadata.filt$Shannon, metadata.filt$annot, kw = TRUE, label = TRUE, method = "bh")
  print(dunn_test)
}

# Summarize data for boxplot visualization
summary_df <- metadata.filt %>%
  group_by(sample.title) %>%
  summarise(
    median_value = median(Shannon),
    Q1 = quantile(Shannon, 0.25),  # First quartile (Q1)
    Q3 = quantile(Shannon, 0.75)   # Third quartile (Q3)
  )

# ----------------------------------------
# Plotting box plot
# ----------------------------------------

# Generate alpha diversity boxplot with significance annotations
alpha.plot = ggplot(metadata.filt, aes(x = annot, y = Shannon, fill = annot)) + 
  geom_boxplot(outliers = FALSE, width = 0.5, position = position_dodge(width = 0.5), color = "black") +
  geom_point(position = position_jitter(width = 0.2), colour = "black", size = 2) +
  theme_prism() +
  geom_signif(comparisons = list(comp1=c("D2 HMA", "WT"),
                                 comp3=c("D7 HMA", "WT"),
                                 comp2=c("D2 HMA", "D7 HMA")),
              map_signif_level = FALSE, step_increase = 0.1, textsize = 4, y_position = c(4.4),
              color = "black") +
  scale_fill_manual(values = c("WT" = "gray43", "D2 HMA" = "steelblue3", "D7 HMA" = "indianred3")) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 2, color = "black") +  # Median as a horizontal line
  theme(axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),  # Increase y-axis label size
        legend.position = "none",
        axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(3.1,4.7))

# Display plot
alpha.plot



