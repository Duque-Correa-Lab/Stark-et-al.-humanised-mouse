# Load required libraries
library(vegan)         # For diversity analysis
library(ggplot2)       # For plotting
library(RColorBrewer)  # For color palettes
library(stringr)       # For string manipulation
library(gridExtra)     # For arranging ggplots
library(dunn.test)     # For Dunn's test (non-parametric post-hoc)
library(ggprism)       # For prism-style ggplots
library(ggsignif)      # For adding significance to plots
library(multcompView)  # For multiple comparisons visualization
library(tidyverse)     # For tidy data manipulation

# ----------------------------------------
# Load and preprocess data
# ----------------------------------------

# Set working directory and load OTU data
cat("Loading data ...\n")
setwd("/your/path/here/") # <-- Use your actual path
data_otu = read.delim("merged_abundance_table.txt", skip=1L)

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
species.df = t(species.df[,-1]) # Transpose and remove clade_name column

# Load and align metadata with species data
metadata = read.delim("AerobevsAnaerobe_scraping_metadata.txt")
rownames(metadata) = metadata$run.accession

# ----------------------------------------
# Alpha diversity analysis (Shannon-index)
# ----------------------------------------

# Calculate Shannon diversity and merge with metadata
cat("Calculating Shannon-index ...\n")
alpha.div = diversity(species.df, index = "shannon")
metadata.filt = metadata[rownames(species.df), ]
metadata.filt$Shannon = alpha.div[rownames(metadata.filt)]
metadata.filt = metadata.filt[metadata.filt$annot1 != "WT", ]
metadata.filt = metadata.filt[metadata.filt$annot1 != "D2 HMA", ]

metadata.filt$annot2 = factor(metadata.filt$annot2, levels = c("Uncultured", "Aerobic", "Anaerobic"))

# Perform Kruskal-Wallis test
kruskal_test_result <- kruskal.test(Shannon ~ annot2, data = metadata.filt)
kruskal_p <- kruskal_test_result$p.value

# Post-hoc Dunn’s test and adjust p-values if Kruskal-Wallis is significant
if (kruskal_p < 0.05) {
  dunn_test <- dunn.test(metadata.filt$Shannon, metadata.filt$annot2, kw = TRUE, label = TRUE)
  dunn_results <- data.frame(dunn_test$res)
  
  # Adjust p-values using Benjamini-Hochberg correction
  dunn_results$p.adj <- p.adjust(dunn_results$P, method = "BH")
  
  # Create a "Comparison" column by concatenating group names
  dunn_results$Comparison <- paste(dunn_results$Group1, dunn_results$Group2, sep = " - ")
  
  # Filter significant comparisons
  significant_pairs <- dunn_results %>%
    filter(p.adj < 0.05) %>%
    select(Comparison, p.adj)
}

# Prepare list of significant pairs for plotting
signif_pairs <- strsplit(as.character(significant_pairs$comparison), " - ")

# Summarize data for boxplot visualization
summary_df <- metadata.filt %>%
  group_by(annot2) %>%
  summarise(
    median_value = median(Shannon),
    Q1 = quantile(Shannon, 0.25),  # First quartile (Q1)
    Q3 = quantile(Shannon, 0.75)   # Third quartile (Q3)
  )

# ----------------------------------------
# Plotting box plot
# ----------------------------------------

# Generate alpha diversity boxplot with significance annotations
alpha.plot = ggplot(metadata.filt, aes(annot2, Shannon, colour = annot1, shape = annot2, fill = annot1)) + 
  geom_boxplot(outliers = FALSE, coef = 0, width = 0.5, position = position_dodge(width = 0.5)) +
  theme_prism() +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  stat_signif(comparisons = signif_pairs,
              map_signif_level = TRUE,
              step_increase = 0.1) +
  labs(x = "", y = "Shannon-index") +
  theme(legend.title = element_blank()) +
  theme(axis.title = element_text(size=14, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold")) +  # Increase y-axis label size
  theme(axis.text.x = element_text(size = 18, face = "bold"),  legend.position = "none") +
  geom_signif(comparisons = list(comp1=c("Aerobic", "Anaerobic"),
                                 comp3=c("Aerobic", "Uncultured"),
                                 comp2=c("Uncultured", "Anaerobic")),
              map_signif_level = FALSE, step_increase = 0.1, textsize = 4, y_position = c(3.2),
              color = "black") +
  scale_colour_manual(values = c("D7 HMA" = "indianred4"))+
  scale_fill_manual(values = c( "D7 HMA" = "indianred4"))+
  scale_shape_manual(values = c(15, 16, 17))

# Display plot
alpha.plot

