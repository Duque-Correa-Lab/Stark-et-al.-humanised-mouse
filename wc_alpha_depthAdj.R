#
# Whole-caecum comparison — WT, D2 HMA, D7 HMA
#
# ============================
# 1. Load libraries
# ============================
library(tidyverse)
library(vegan)
library(ggprism)      
library(ggsignif)     
library(dunn.test)    

setwd("your path")

# Standard asterisk convention for p-values, shared across plots in this script

p_to_stars <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.0001) "****"
  else if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

# ============================
# 2. Read files
# ============================
meta_wc <- read_tsv("metadata_WC_HMAvsWT.txt")
seq_wc  <- read_csv("read_counts_whole_caecum.csv")
otu_wc  <- read.delim("merged_abundance_table_whole_caecum.txt", skip = 1, check.names = FALSE)

# ============================
# 3. Keep species only
# ============================
otu_wc_sp <- otu_wc %>%
  filter(grepl("s__", clade_name), !grepl("t__", clade_name)) %>%
  distinct(clade_name, .keep_all = TRUE)

# Set taxa as rownames
otu_mat_wc <- otu_wc_sp %>%
  column_to_rownames("clade_name") %>%
  as.matrix()

# ============================
# 4. Transpose: samples × species
# ============================
otu_mat_wc <- t(otu_mat_wc)

# Remove _results_gtdb suffix from sample names
rownames(otu_mat_wc) <- sub("_results_gtdb", "", rownames(otu_mat_wc))

# ============================
# 5. Match sample order between read-count file and OTU table
# ============================
# Keep only samples present in OTU table
seq_wc_filt <- seq_wc %>%
  filter(sample %in% rownames(otu_mat_wc))

# Reorder to match OTU matrix
seq_wc_filt <- seq_wc_filt[match(rownames(otu_mat_wc), seq_wc_filt$sample), ]

# Confirm alignment
stopifnot(all(rownames(otu_mat_wc) == seq_wc_filt$sample))

# ============================
# 6. Reconstruct pseudo-counts from relative abundance x total reads
# ============================
otu_mat_prop <- otu_mat_wc / 100
otu_mat_counts <- otu_mat_prop
for(i in 1:nrow(otu_mat_counts)){
  otu_mat_counts[i, ] <- round(otu_mat_prop[i, ] * seq_wc_filt$reads_1[i])
}
otu_mat_counts[is.na(otu_mat_counts)] <- 0
mode(otu_mat_counts) <- "integer"

# ============================
# 7. Rarefy to the minimum reconstructed species-level read total
# ============================

otu_mat_counts_nonzero <- otu_mat_counts[rowSums(otu_mat_counts) > 0, , drop = FALSE]
min_depth <- min(rowSums(otu_mat_counts_nonzero))
set.seed(123)
otu_mat_rarefied_wc <- vegan::rrarefy(otu_mat_counts_nonzero, sample = min_depth)

# ============================
# 8. Shannon diversity
# ============================
shannon_wc <- diversity(otu_mat_rarefied_wc, index = "shannon")
shannon_df_wc <- tibble(
  sample  = rownames(otu_mat_rarefied_wc),
  Shannon = shannon_wc
)

# ============================
# 9. Merge metadata
# ============================
shannon_df_wc_clean <- shannon_df_wc %>%
  rename_with(~ "run_accession", "sample")

shannon_meta_wc <- shannon_df_wc_clean %>%
  left_join(meta_wc, by = "run_accession")

# View result
head(shannon_meta_wc)

# ============================
# 10. Kruskal-Wallis test
# ============================
kruskal_test_result <- kruskal.test(Shannon ~ annot, data = shannon_meta_wc)
kruskal_test_result

# If significant, run Dunn's test and extract the WT vs D7 HMA comparison;
# otherwise, skip Dunn's test and label the plot as non-significant.

comparison_name <- "D7 HMA - WT"
show_signif <- FALSE

if (kruskal_test_result$p.value < 0.05) {
  dunn_test <- dunn.test(shannon_meta_wc$Shannon, shannon_meta_wc$annot, kw = TRUE, label = TRUE, method = "bh")
  print(dunn_test)

  p_value <- dunn_test$P.adjusted[
  which(dunn_test$comparisons == comparison_name)
  ]

  # Format label
  p_label <- paste0("p = ", signif(p_value, 3))
  p_label

  p_label <- ifelse(p_value < 0.0001, "***",
                  ifelse(p_value < 0.001, "**",
                         ifelse(p_value < 0.05, "*", "ns")))
  show_signif <- p_label != "ns"
  
} else {
  
  message("Kruskal-Wallis test not significant (p = ", signif(kruskal_test_result$p.value, 3),
          "). Skipping Dunn's test and significance bracket.")
  
  p_value <- NA
  p_label <- "ns"
  
}

# ============================
# 11. Summary statistics
# ============================
summary_df <- shannon_meta_wc %>%
  group_by(annot) %>%
  summarise(
    mean_value = mean(Shannon),
    sd_value   = sd(Shannon)
  )
summary_df

# ============================
# 12. Plot: Boxplot + median crossbar + jittered points
# ============================
# Set factor order for x-axis
shannon_meta_wc <- shannon_meta_wc %>%
  mutate(annot = factor(annot, levels = c("WT", "D2 HMA", "D7 HMA")))

# Define colors
colors <- c("WT" = "gray43", "D2 HMA" = "steelblue3", "D7 HMA" = "indianred3")

#Significance bar position relative to WT only
wt_max <- max(shannon_meta_wc$Shannon[shannon_meta_wc$annot == "WT"])

# Plot

alpha.plot <- ggplot(shannon_meta_wc,
                     aes(x = annot, y = Shannon, fill = annot, color = annot)) +
  
  # Layer 1: box outline only (black)
  geom_boxplot(
    color = "black",
    width = 0.5,
    outlier.shape = NA,
    alpha = 0.7
  ) +

  # Layer 2: median line colored by group
  stat_summary(
    aes(color = annot),
    fun = median,
    fun.min = median,
    fun.max = median,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.8
  ) +
  
  # Points
  geom_point(
    aes(color = annot),
    position = position_jitter(width = 0.25),
    size = 6,
    show.legend = FALSE
  ) +
  
  # Manual colors
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  
  theme_prism() +
  
  labs(x = "", y = "Shannon-index") +
  
  theme(
    axis.title = element_text(size = 26, face = "bold"),
    axis.text  = element_text(size = 26, face = "bold"),
    axis.title.y = element_text(size = 32, face = "bold"),
    legend.position = "none"
  )
  
# Significance annotation using stars (don't show NS)
if (show_signif) {
  alpha.plot <- alpha.plot +
    geom_signif(
    comparisons = list(c("WT", "D7 HMA")),
    annotations = p_label,
    y_position = wt_max + 0.05,
    tip_length = 0,
    size = 1.2,
    textsize = 16,
    color = "black"
  ) 
}

# ============================
# 13. Show plot
# ============================

alpha.plot

ggsave(alpha.plot, filename="alpha-div_WC_depthAdj.pdf", height=14, width=14, dpi=300)
