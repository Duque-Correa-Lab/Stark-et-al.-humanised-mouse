# =============================================================================
# BLOCK 1: Load libraries
# =============================================================================

library(file2meco)   # provides humann2meco() and MetaCyc_pathway_map
library(microeco)    # provides microtable and trans_func
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(Maaslin2)
library(RColorBrewer)
library(pheatmap)
library(grid)
library(ggpubr)
library(ggprism)


# =============================================================================
# BLOCK 2: Load input data
# =============================================================================

# Load the MetaCyc hierarchy map (Superclass1, Superclass2, Pathway columns)
data("MetaCyc_pathway_map")
head(MetaCyc_pathway_map)
colnames(MetaCyc_pathway_map)

setwd("/Users/kas206/Library/Mobile Documents/com~apple~CloudDocs/Documents/Paper_revision/Metagenomics/GitHub/") # edit with real path

# Load sample metadata 
sample_meta <- read.delim("metadata_WC_SC_C.txt", row.names = 1)

# Convert HUMAnN merged pathway abundance table to a microtable object,
# incorporating sample metadata
tmp_df <- read.delim("pathabundances_merged.tsv", check.names = FALSE)
colnames(tmp_df) <- gsub("_merged_Abundance", "", colnames(tmp_df))
any(is.na(tmp_df))
tmp_df[is.na(tmp_df)] <- 0
write.table(tmp_df, "pathabundances_cleaned.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
mt <- humann2meco(
  feature_table = "pathabundances_cleaned.tsv",
  sample_table  = sample_meta,
  db            = "MetaCyc"
)

# Inspect the microtable
mt
mt$tax_table |> head()   # Superclass1, Superclass2, Pathway columns
mt$otu_table |> head()   # pathway abundances (pathways x samples)


# =============================================================================
# BLOCK 3: Clean the microtable
# =============================================================================

# Strip suffix added by HUMAnN during sample merging
colnames(mt$otu_table) <- gsub("_merged_Abundance", "", colnames(mt$otu_table))

# Restrict to the D7 HMA timepoint, luminal content vs mucosal scraping only.
# mt$sample_table$annot1 has 3 levels (WT, D2 HMA, D7 HMA) and annot2 has 3
# levels (luminal content, mucosal scraping, whole caecum) - this script is
# specifically the D7 HMA scraping-vs-content comparison, so both 
# filtered explicitly rather than relying on factor() to do it (factor()
# only affects display/level order - it does NOT remove non-matching rows,
# it silently turns them into NA while leaving them in the dataset).

mt$sample_table <- mt$sample_table[
  mt$sample_table$annot1 == "D7 HMA" &
  mt$sample_table$annot2 %in% c("luminal content", "mucosal scraping"),
]
mt$tidy_dataset()   # sync otu_table/tax_table to the now-restricted sample_table

# Keep sample_meta (used directly as MaAsLin2 input_metadata below) in sync
# with the same restriction applied to mt$sample_table
sample_meta <- sample_meta[
  sample_meta$annot1 == "D7 HMA" &
  sample_meta$annot2 %in% c("luminal content", "mucosal scraping"),
]

# Set group factor levels to control display order throughout all figures
mt$sample_table$annot2 <- factor(
  mt$sample_table$annot2,
  levels = c("luminal content", "mucosal scraping")
)

# Remove unclassified, unmapped and ungrouped pathways from the taxonomy table,
# then sync the removal through the abundance and sample tables
mt$tax_table <- mt$tax_table[
  !grepl("unclassified|UNMAPPED|UNGROUPED", mt$tax_table$pathway, ignore.case = TRUE),
]
mt$tidy_dataset()


# =============================================================================
# BLOCK 4: Calculate and export abundance tables
# =============================================================================

# Pre-compute relative abundances at all taxonomic levels (Superclass1,
# Superclass2, pathway). Results stored in mt$taxa_abund.
mt$cal_abund()

# Extract each level into its own data frame (rows = features, columns = samples)
pathway_abund    <- mt$taxa_abund$pathway
superclass1_abund <- mt$taxa_abund$Superclass1
superclass2_abund <- mt$taxa_abund$Superclass2

# Save abundance tables for supplementary materials.
# cbind() prepends rownames as an explicit named column before writing.
write.table(
  cbind(pathway = rownames(pathway_abund), pathway_abund),
  file      = "pathway_abundance_D7HMA.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

write.table(
  cbind(Superclass1 = rownames(superclass1_abund), superclass1_abund),
  file      = "metacyc_superclass1_abundance_D7HMA.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

write.table(
  cbind(Superclass2 = rownames(superclass2_abund), superclass2_abund),
  file      = "metacyc_superclass2_abundance_D7HMA.tsv",
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# Also save long-format versions (one row per feature-sample combination)
# for downstream plotting and analysis
superclass1_long <- superclass1_abund |>
  rownames_to_column("Superclass1") |>
  pivot_longer(-Superclass1, names_to = "Sample", values_to = "Abundance") |>
  arrange(Sample, desc(Abundance))

write.csv(superclass1_long, "metacyc_superclass1_long_D7HMA.csv", row.names = FALSE)

superclass2_long <- superclass2_abund |>
  rownames_to_column("Superclass2") |>
  pivot_longer(-Superclass2, names_to = "Sample", values_to = "Abundance") |>
  arrange(Sample, desc(Abundance))

write.csv(superclass2_long, "metacyc_superclass2_long_D7HMA.csv", row.names = FALSE)


# =============================================================================
# BLOCK 5: Bar charts — top 15 Superclass2 and pathway (all samples)
# =============================================================================

# Build trans_abund objects for Superclass2 and pathway level,
# averaging within each group defined by annot2
t1 <- trans_abund$new(
  dataset   = mt,
  taxrank   = "Superclass2",
  ntaxa     = 15,
  groupmean = "annot2"
)

t2 <- trans_abund$new(
  dataset   = mt,
  taxrank   = "pathway",
  ntaxa     = 15,
  groupmean = "annot2"
)

my_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(16)

# Top 15 Superclass2 bar chart
p_super2 <- t1$plot_bar(
  others_color       = "grey80",
  legend_text_italic = FALSE,
  color_values = my_colors
) +
  geom_bar(stat = "identity", width = 0.6) +
  ylab("Relative Abundance of Classified Pathways (%)") +
  guides(fill = guide_legend(title = "Superpathway")) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggsave("pathway_abundance_D7HMA_superclass2_cleaned.pdf", plot = p_super2, width = 10, height = 6)
ggsave("pathway_abundance_D7HMA_superclass2_cleaned.png", plot = p_super2, width = 10, height = 6, dpi = 300)

# Top 15 pathway bar chart
p_path <- t2$plot_bar(
  others_color       = "grey70",
  legend_text_italic = FALSE,
  color_values = my_colors
) +
  ylab("Relative Abundance of Classified Pathways (%)") +
  ggtitle("Top 15 Classified MetaCyc Pathways by annot2")

ggsave("pathway_abundance_D7HMA_pathway_cleaned.pdf", plot = p_path, width = 8, height = 6)
ggsave("pathway_abundance_D7HMA_pathway_cleaned.png", plot = p_path, width = 8, height = 6, dpi = 300)


# =============================================================================
# BLOCK 6: MaAsLin2 — differential abundance analysis
# All three levels tested with TSS normalisation and LOG transformation.
# Comparison is luminal content vs mucosal scraping (annot2) within the
# D7 HMA timepoint (mt was restricted to D7 HMA in Block 3). Reference
# group is "mucosal scraping", so fold changes are expressed relative to
# mucosal scraping.
# =============================================================================

fit_superclass1 <- Maaslin2(
  input_data     = superclass1_abund,
  input_metadata = sample_meta,
  output         = "maaslin_superclass1_D7HMA",
  normalization  = "TSS",
  transform      = "LOG",
  fixed_effects  = c("annot2"),
  random_effects = NULL,
  reference      = c("annot2", "mucosal scraping"),
  min_prevalence = 0.1,
  max_significance = 0.05
)

fit_superclass2 <- Maaslin2(
  input_data     = superclass2_abund,
  input_metadata = sample_meta,
  output         = "maaslin_superclass2_D7HMA",
  normalization  = "TSS",
  transform      = "LOG",
  fixed_effects  = c("annot2"),
  random_effects = NULL,
  reference      = c("annot2", "mucosal scraping"),
  min_prevalence = 0.1,
  max_significance = 0.05
)

fit_pathway <- Maaslin2(
  input_data     = pathway_abund,
  input_metadata = sample_meta,
  output         = "maaslin_pathway_D7HMA",
  normalization  = "TSS",
  transform      = "LOG",
  fixed_effects  = c("annot2"),
  random_effects = NULL,
  reference      = c("annot2", "mucosal scraping"),
  min_prevalence = 0.1,
  max_significance = 0.05
)


# =============================================================================
# BLOCK 7: Pathways unique to mucosal scraping (Venn analysis)
# =============================================================================

# Merge samples by group, then identify pathways present only in mucosal scraping
# but absent from both WT and D2 HMA
group_mt <- mt$merge_samples(group = "annot2")
t3       <- trans_venn$new(dataset = group_mt)
venn_data <- t3$venn_list

target         <- "mucosal scraping"
others         <- names(venn_data)[names(venn_data) != target]
unique_to_mucosa <- setdiff(venn_data[[target]], unlist(venn_data[others]))
print(unique_to_mucosa)


# =============================================================================
# BLOCK 8: Extract significant features from MaAsLin2 results
# MaAsLin2 (via R's make.names()) encodes feature names by converting any
# non-alphanumeric character to a dot, and prepends "X" to names that don't
# start with a letter. Both the feature names and the abundance table
# rownames are normalized to a common format (see clean_feature_name() /
# strip_artifact_X() below) before matching.
# =============================================================================

# make.names()/MaAsLin prepends "X" to feature names that don't start with
# a letter (e.g. names beginning with a digit, or a symbol that gets
# sanitized to a dot). Strip that artifact - but only when "X" is followed
# by a digit, a dot, or a space, never by another letter, so genuine
# pathway names that legitimately start with "X" (e.g. "Xylose...") are
# left alone.
strip_artifact_X <- function(x) {
  gsub("^X(?=[0-9]|\\.| )", "", x, perl = TRUE)
}

# Shared cleaning function applied identically to MaAsLin feature names and
# pathway_abund rownames. Rather than whitelist individual punctuation
# characters one at a time, collapse ANY run of non-alphanumeric characters
# into a single space - this mirrors what MaAsLin's own sanitization already
# does to build feature names, so both sides end up in the same normalized
# form regardless of which specific punctuation/HTML entity/Greek letter a
# given pathway name happens to contain.

clean_feature_name <- function(x) {
  x <- gsub("\\&\\&", " AND ", x)      # preserve the AND relationship in "&&"-joined superpathway names before it gets stripped below
  x <- gsub("[^A-Za-z0-9]+", " ", x)   # collapse any run of non-alphanumeric characters to one space
  trimws(x)
}

# --- Significant pathways ---
res_pathway     <- read.delim("maaslin_pathway_D7HMA/all_results.tsv")
sig_pathway     <- res_pathway |> dplyr::filter(metadata == "annot2", qval < 0.05, abs(coef) > 2 )

sig_pathway_clean <- sig_pathway$feature |> strip_artifact_X() |> clean_feature_name()

pathway_abund_rownames <- clean_feature_name(rownames(pathway_abund))

intersect(sig_pathway_clean, pathway_abund_rownames)   # verify matches

matched_idx <- pathway_abund_rownames %in% sig_pathway_clean
df_sig_p <- pathway_abund[matched_idx, ] |>
  as.data.frame() |>
  rownames_to_column("pathway_clean") |>
  pivot_longer(-pathway_clean, names_to = "Sample", values_to = "Abundance") |>
  left_join(mt$sample_table |> rownames_to_column("Sample"), by = "Sample")

# --- Significant Superclass1 features ---
# NOTE: this claim was written for the WT/D2 HMA/D7 HMA whole-caecum
# comparison and has not been reverified for this D7-HMA scraping-vs-content
# comparison - check nrow(sig_superclass1) for this dataset before relying
# on it. df_sig_s1 is retained either way as a template for downstream use.

res_superclass1  <- read.delim("maaslin_superclass1_D7HMA/all_results.tsv")
sig_superclass1  <- res_superclass1 |> dplyr::filter(metadata == "annot2", qval < 0.05)

sig_superclass1_clean <- sig_superclass1$feature |> strip_artifact_X() |> clean_feature_name()

superclass1_abund_rownames <- clean_feature_name(rownames(superclass1_abund))

intersect(sig_superclass1_clean, superclass1_abund_rownames)   # verify matches

matched_idx <- superclass1_abund_rownames %in% sig_superclass1_clean
df_sig_s1 <- superclass1_abund[matched_idx, ] |>
  as.data.frame() |>
  rownames_to_column("Superclass1") |>
  pivot_longer(-Superclass1, names_to = "Sample", values_to = "Abundance") |>
  left_join(mt$sample_table |> rownames_to_column("Sample"), by = "Sample")

# --- Significant Superclass2 features ---
res_superclass2  <- read.delim("maaslin_superclass2_D7HMA/all_results.tsv")
sig_superclass2  <- res_superclass2 |> dplyr::filter(metadata == "annot2", qval < 0.05)

sig_superclass2_clean <- sig_superclass2$feature |> strip_artifact_X() |> clean_feature_name()

superclass2_abund_rownames <- clean_feature_name(rownames(superclass2_abund))

intersect(sig_superclass2_clean, superclass2_abund_rownames)   # verify matches

matched_idx <- superclass2_abund_rownames %in% sig_superclass2_clean
df_sig_s2 <- superclass2_abund[matched_idx, ] |>
  as.data.frame() |>
  rownames_to_column("Superclass2") |>
  pivot_longer(-Superclass2, names_to = "Sample", values_to = "Abundance") |>
  left_join(mt$sample_table |> rownames_to_column("Sample"), by = "Sample")
df_sig_s2 <- df_sig_s2 %>%
  mutate(Abundance = as.numeric(Abundance))

# =============================================================================
# BLOCK 9: Bar chart — significant Superclass2 features
# Abundances from cal_abund() are already relative (0-1 scale), so multiply
# by 100 to express as percentages. Bars therefore show the true percentage
# of classified reads, not a re-normalised fraction of significant terms only.
# =============================================================================

df_sig_rel_s2 <- df_sig_s2 |>
  mutate(RelAbundance_new = Abundance * 100)

# Sanity check: total plotted abundance per group (informational only - not
# expected to sum to 100%, since only significant Superclass2 terms are here)
df_sig_rel_s2 %>%
  group_by(annot2) %>%
  summarise(total = sum(RelAbundance_new, na.rm = TRUE))

colors_16   <- brewer.pal(n = 12, name = "Paired")
colors_extra <- brewer.pal(n = 8,  name = "Set2")[1:4]
colors_16   <- c(colors_16, colors_extra)

pbar <- ggplot(df_sig_rel_s2, aes(x = annot2, y = RelAbundance_new, fill = Superclass2)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colors_16) +
  labs(
    y    = "Relative Abundance of Classified Pathways (%)",
    x    = "",
    fill = "Superclass2"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 0, hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black"),
    legend.text  = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold")
  )

pbar

ggsave("superclass2_abundance_D7HMA_NOunclassified_sig.pdf", plot = pbar, width = 8, height = 6)
ggsave("superclass2_abundance_D7HMA_NOunclassified_sig.png", plot = pbar, width = 8, height = 6, dpi = 300)


# =============================================================================
# BLOCK 10: Heatmaps — significant pathway features
# Two options are provided: a ggplot2 tile heatmap (Option 1) and a
# pheatmap clustered heatmap with sample annotations (Option 2).
# Both use df_sig_p built in Block 8.
# =============================================================================

# Restore factor order after joins (ensures correct panel order in facet_grid)
df_sig_p$annot2 <- factor(df_sig_p$annot2, levels = c("luminal content", "mucosal scraping"))

# Order pathways by total abundance for consistent row ordering
pathway_order <- df_sig_p |>
  group_by(pathway_clean) |>
  summarise(total_abundance = sum(Abundance)) |>
  arrange(desc(total_abundance)) |>
  pull(pathway_clean)

df_sig_p$pathway_clean <- factor(df_sig_p$pathway_clean, levels = pathway_order)

# --- Option 1: ggplot2 tile heatmap ---
hmap <- ggplot(df_sig_p, aes(x = Sample, y = pathway_clean, fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  facet_grid(~ annot2, scales = "free_x", space = "free") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y  = element_text(size = 10, face = "bold"),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    strip.text   = element_text(size = 9, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10)
  ) +
  labs(fill = "Abundance")

hmap

ggsave("heatmap_significant_pathways.png", hmap, width = 12, height = 8, dpi = 300)
ggsave("heatmap_significant_pathways.pdf", hmap, width = 12, height = 8)

# --- Option 2: pheatmap with hierarchical clustering and sample annotations ---

# Reshape to pathways x samples matrix
data_otu <- df_sig_p |>
  select(Sample, pathway_clean, Abundance) |>
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0)

data_otu_mat <- as.matrix(data_otu[, -1])
rownames(data_otu_mat) <- data_otu$pathway_clean

# Sample annotation bar colours
sample_coloring <- df_sig_p |>
  select(Sample, annot2) |>
  distinct() |>
  column_to_rownames("Sample")

palette_annot <- c("luminal content" = "#FA8", "mucosal scraping" = "#9E2A2F")
ann_colors    <- list(annot2 = palette_annot)

# Plot pheatmap with row z-score scaling and clustering on both axes
phmap <- pheatmap(
  mat                  = data_otu_mat,
  scale                = "row",
  show_rownames        = TRUE,
  show_colnames        = FALSE,
  cluster_rows         = TRUE,
  cluster_cols         = TRUE,
  annotation_col       = sample_coloring,
  annotation_names_col = FALSE,
  annotation_colors    = ann_colors,
  fontsize_row         = 4,
  cellwidth            = 10,
  cellheight           = 4
)

# Use grid.draw rather than ggsave — pheatmap does not produce a ggplot object
png("heatmap_significant_pathways_pheatmap.png", width = 9, height = 10, units = "in", res = 300)
grid::grid.draw(phmap$gtable)
dev.off()

pdf("heatmap_significant_pathways_pheatmap.pdf", width = 9, height = 10)
grid::grid.draw(phmap$gtable)
dev.off()

# Highlight the chitin deacetylation row label in red
row_labels  <- phmap$gtable$grobs[[which(phmap$gtable$layout$name == "row_names")]]
chitin_idx  <- which(row_labels$label == "chitin deacetylation")
label_colours <- rep("black", length(row_labels$label))
label_colours[chitin_idx] <- "red"
row_labels$gp$col <- label_colours
phmap$gtable$grobs[[which(phmap$gtable$layout$name == "row_names")]] <- row_labels

png("heatmap_significant_pathways_chitin.png", width = 9, height = 10, units = "in", res = 300)
grid::grid.draw(phmap$gtable)
dev.off()

pdf("heatmap_significant_pathways_chitin.pdf", width = 9, height = 10)
grid::grid.draw(phmap$gtable)
dev.off()


# =============================================================================
# BLOCK 11: Single pathway plot — "chitin deacetylation"
# Individual sample jitter with mean and SEM.
# =============================================================================

grep("chitin deacetylation", rownames(pathway_abund), value = TRUE)

pathway_name <- "chitin deacetylation"

df_path <- pathway_abund[pathway_name, , drop = FALSE] |>
  as.data.frame() |>
  rownames_to_column("pathway") |>
  pivot_longer(-pathway, names_to = "Sample", values_to = "Abundance")


sample_meta_fixed <- sample_meta |>
  rownames_to_column("Sample")

df_path <- df_path |>
  left_join(sample_meta_fixed, by = "Sample")


# --- MaAsLin2 Significance Extraction ---
q_to_stars <- function(q) {
  if (length(q) == 0 || is.na(q)) return("ns")
  if (q < 0.0001) "***"
  else if (q < 0.001) "**"
  else if (q < 0.05) "*"
  else "ns"
}

# Match naming variations in your res_pathway tracking
chitin_row <- res_pathway |>
  dplyr::filter(
    (feature == "chitin.deacetylation"), 
    metadata == "annot2"
  )

# Extract q-value and map to significance star string
q_val_maaslin     <- chitin_row$qval
sig_label_maaslin <- q_to_stars(q_val_maaslin)

# Dynamic bracket sizing calculations
ymax   <- max(df_path$Abundance, na.rm = TRUE)
yrange <- diff(range(df_path$Abundance, na.rm = TRUE))
if (yrange == 0) yrange <- ymax * 0.1 


# Plot

plot_chitin <- ggplot(df_path, aes(x = annot2, y = Abundance, color = annot2)) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", size = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, color = "black", size = 0.8) +
  scale_y_continuous(
    limits  = c(0, NA),
    expand  = expansion(mult = c(0, 0.20)) # Headroom for the significance bar
  ) +
  scale_color_manual(values = c("luminal content" = "#FA8", "mucosal scraping" = "#9E2A2F")) +
  # Replaced stat_compare_means with geom_signif using MaAsLin2 metrics
  ggsignif::geom_signif(
    comparisons = list(c("luminal content", "mucosal scraping")),
    annotations = sig_label_maaslin,
    y_position  = ymax + yrange * 0.08,
    tip_length  = 0.01,
    bracket.size = 0.8,
    size        = 1.2,
    textsize    = 6,
    color       = "black"
  ) +
  labs(
    title = "Chitin deacetylation pathway",
    x     = "",
    y     = "Abundance"
  ) +
  theme_prism(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text       = element_text(size = 12, face = "bold"),
    axis.title      = element_text(size = 14, face = "bold")
  )

plot_chitin

ggsave("chitin_deacetylation_D7HMA.pdf", plot = plot_chitin, width = 5, height = 5)
ggsave("chitin_deacetylation_D7HMA.png", plot = plot_chitin, width = 5, height = 5, dpi = 300)


# =============================================================================
# BLOCK 12: Bar chart — absolute pathway abundance by group
# Unlike the relative abundance plots above, this uses raw otu_table values
# aggregated by Superclass2, showing mean absolute abundance per group.
# =============================================================================

# Build long-format data from otu_table and join taxonomy and metadata
df_long <- mt$otu_table |>
  as.data.frame() |>
  rownames_to_column("pathway") |>
  pivot_longer(-pathway, names_to = "Sample", values_to = "Abundance")

tax_df <- mt$tax_table |>
  rownames_to_column("pathway_id") |>
  rename(pathway_name = pathway)

df_long <- df_long |>
  left_join(tax_df, by = c("pathway" = "pathway_id")) |>
  left_join(mt$sample_table |> rownames_to_column("Sample"), by = "Sample")

# Aggregate to Superclass2 level, collapse all but top 15 into "Others"
df_grouped <- df_long |>
  group_by(Sample, annot2, Superclass2) |>
  summarise(Abundance = sum(Abundance), .groups = "drop")

top15 <- df_grouped |>
  group_by(Superclass2) |>
  summarise(Total = sum(Abundance)) |>
  arrange(desc(Total)) |>
  slice(1:15) |>
  pull(Superclass2)

df_grouped <- df_grouped |>
  mutate(Superclass2 = ifelse(Superclass2 %in% top15, Superclass2, "Others")) |>
  group_by(Sample, annot2, Superclass2) |>
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Average across samples within each group
df_mean <- df_grouped |>
  group_by(annot2, Superclass2) |>
  summarise(Abundance = mean(Abundance), .groups = "drop")

df_mean$Superclass2 <- factor(df_mean$Superclass2, levels = c(top15, "Others"))

p_abs <- ggplot(df_mean, aes(x = annot2, y = Abundance, fill = Superclass2)) +
  geom_bar(stat = "identity") +
  ylab("Mean Absolute Pathway Abundance") +
  xlab("") +
  ggtitle("MetaCyc Pathways (Absolute Abundance)") +
  theme_bw() +
  theme(
    axis.text.x  = element_text(size = 12),
    legend.title = element_blank()
  )

ggsave("pathway_abundance_D7HMA_absolute.pdf", plot = p_abs, width = 7, height = 5)
ggsave("pathway_abundance_D7HMA_absolute.png", plot = p_abs, width = 7, height = 5, dpi = 300)
