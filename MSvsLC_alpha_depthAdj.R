# Compares Shannon index for mucosal vs luminal samples (from WT, D2 HMA and D7 HMA)
#
# ============================
# 1. Load libraries
# ============================
library(tidyverse)
library(vegan)
library(ggprism)      # for theme_prism()
library(ggsignif)     # for geom_signif()
library(dunn.test)    # for Dunn post-hoc
library(dplyr)
library(ggplot2)

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

meta_scraping <- read_tsv("AerobevsAnaerobe_scraping_metadata.txt")
seq_scraping  <- read_csv("read_counts_scraping.csv")
otu_scraping  <- read.delim("merged_abundance_table_scraping.txt", skip = 1, check.names = FALSE)

meta_content <- read_tsv("AerobevsAnaerobe_content_metadata.txt")
seq_content  <- read_csv("read_counts_content.csv")
otu_content  <- read.delim("merged_abundance_table_content.txt", skip = 1, check.names = FALSE)

prepare_otu <- function(otu_raw) {
  otu_sp <- otu_raw %>%
    filter(grepl("s__", clade_name), !grepl("t__", clade_name)) %>%
    distinct(clade_name, .keep_all = TRUE)
  
  otu_mat <- otu_sp %>%
    column_to_rownames("clade_name") %>%
    as.matrix()
  
  otu_mat <- t(otu_mat)
  rownames(otu_mat) <- sub("_results_gtdb", "", rownames(otu_mat))
  
  return(otu_mat)
}

otu_scraping_mat <- prepare_otu(otu_scraping)
otu_content_mat  <- prepare_otu(otu_content)



# ============================
# 3. Match sample order between read-count files and OTU tables
# ============================

align_depth <- function(otu_mat, seq_df) {
  # CSV uses "Run-Accession" (matches OTU matrix
  # rownames, which are run accessions from prepare_otu()) and "Read-count"
  # for total reads - different to the whole caecum script's "sample"/
  # "reads_1" naming. Standardizing to "sample"/"reads_1" internally so the
  # rest of this pipeline (convert_to_counts, etc.) doesn't need to change.
  
  seq_df <- seq_df %>%
    rename(sample = `Run-Accession`, reads_1 = `Read-count`)
  
  shared <- intersect(rownames(otu_mat), seq_df$sample)
  
  seq_filt <- seq_df %>%
    filter(sample %in% shared) %>%
    arrange(match(sample, shared))
  
  otu_mat <- otu_mat[seq_filt$sample, ]
  
  stopifnot(identical(rownames(otu_mat), seq_filt$sample))
  
  return(list(otu = otu_mat, seq = seq_filt))
}

scraping_aligned <- align_depth(otu_scraping_mat, seq_scraping)
content_aligned  <- align_depth(otu_content_mat, seq_content)

otu_scraping_mat <- scraping_aligned$otu
seq_scraping_filt <- scraping_aligned$seq

otu_content_mat <- content_aligned$otu
seq_content_filt <- content_aligned$seq

# ============================
# 4. Convert relative abundance → counts
# ============================

convert_to_counts <- function(otu_mat, seq_df) {
  otu_prop <- otu_mat / 100
  
  otu_counts <- otu_prop
  for(i in 1:nrow(otu_counts)) {
    otu_counts[i, ] <- round(otu_prop[i, ] * seq_df$reads_1[i])
  }
  
  otu_counts[is.na(otu_counts)] <- 0
  mode(otu_counts) <- "integer"
  
  return(otu_counts)
}

otu_scraping_counts <- convert_to_counts(otu_scraping_mat, seq_scraping_filt)
otu_content_counts  <- convert_to_counts(otu_content_mat, seq_content_filt)


# ============================
# 5. Merge OTU tables + metadata
# ============================

scraping_df <- as.data.frame(otu_scraping_counts) %>%
  rownames_to_column("run_accession")

content_df <- as.data.frame(otu_content_counts) %>%
  rownames_to_column("run_accession")

otu_all <- bind_rows(scraping_df, content_df)

meta_scraping$source <- "scraping"
meta_content$source <- "content"

meta_all <- bind_rows(meta_scraping, meta_content)

full_data <- meta_all %>%
  inner_join(otu_all, by = "run_accession")


# ============================
# 6. Filter uncultured samples
# ============================

uncultured <- full_data %>%
  filter(tolower(annot2) == "uncultured")


# ============================
# 7. Rarefy within annot1
# ============================
otu_cols <- colnames(otu_all)[-1]

groups <- uncultured %>% group_split(annot1)

rarefied_list <- lapply(groups, function(df) {
  
  otu <- df[, otu_cols] %>% as.matrix()
  
  otu[is.na(otu)] <- 0
  mode(otu) <- "integer"
  
  # Remove zero-depth samples
  depth <- rowSums(otu)
  otu <- otu[depth > 0, , drop = FALSE]
  df  <- df[depth > 0, ]
  
  # Skip groups with too few samples
  if(nrow(otu) < 2) return(NULL)
  
  # Rarefy
  min_depth <- min(rowSums(otu))
  
  rare <- rrarefy(otu, sample = min_depth) %>% as.data.frame()
  
  # Add metadata back
  rare$run_accession <- df$run_accession
  rare$annot1 <- df$annot1
  rare$source <- df$source   # content vs scraping (if present)
  
  return(rare)
})

otu_rarefied <- bind_rows(rarefied_list)


# ============================
# 8. Compute Shannon
# ============================
shannon <- diversity(otu_rarefied[, otu_cols], index = "shannon")

shannon_meta <- tibble(
  run_accession = otu_rarefied$run_accession,
  annot1 = otu_rarefied$annot1,
  source = otu_rarefied$source,
  Shannon = shannon
)


# ============================
# 9. Plot Shannon MS vs LC within each annot1 
# ============================

plot_configs <- list(
  
  list(
    group_val = "WT",
    colors = c("scraping" = "gray12",
               "content" = "gray50"),
    filename = "Shannon_WT_MSvsLC.pdf",
    title = "WT: Mucosal vs Luminal"
  ),
  
  list(
    group_val = "D2 HMA",
    colors = c("scraping" = "steelblue4",
               "content" = "steelblue2"),
    filename = "Shannon_D2HMA_MSvsLC.pdf",
    title = "D2 HMA: Mucosal vs Luminal"
  ),
  
  list(
    group_val = "D7 HMA",
    colors = c("scraping" = "indianred4",
               "content" = "#FA8"),
    filename = "Shannon_D7HMA_MSvsLC.pdf",
    title = "D7 HMA: Mucosal vs Luminal"
  )
)

# ============================
# Shared y-axis range across all panels (so plots are directly comparable
# side by side for publication)
# ============================
overall_ymin <- min(shannon_meta$Shannon, na.rm = TRUE)
overall_ymax <- max(shannon_meta$Shannon, na.rm = TRUE)
overall_yrange <- overall_ymax - overall_ymin

# Headroom reserved for the single significance bracket each panel can have
y_axis_limits <- c(
  overall_ymin - overall_yrange * 0.05,
  overall_ymax + overall_yrange * 0.2
)

for(cfg in plot_configs) {
  
  df <- shannon_meta %>%
    filter(annot1 == cfg$group_val) %>%
    mutate(source = factor(source,
                           levels = c("content", "scraping")))
  
  # ============================
  # Wilcoxon test
  # ============================
  wt <- wilcox.test(Shannon ~ source, data = df)
  print(cfg$title)
  print(wt)
  
  p_value <- wt$p.value
  
  p_label <- p_to_stars(p_value)


  # ============================
  # Plot
  # ============================

  cont_max <- shannon_meta %>%
    dplyr::filter(
      annot1 == cfg$group_val,
      source == "content"
    ) %>%
    dplyr::pull(Shannon) %>%
    max(na.rm = TRUE)
  
  p <- ggplot(df,
              aes(x = source,
                  y = Shannon,
                  fill = source,
                  color = source)) +
    
    geom_boxplot(
      color = "black",
      width = 0.55,
      alpha = 0.8,
      outlier.shape = NA
    ) +
    
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      color = "black",
      width = 0.55,
      linewidth = 0.6
    ) +
    
    geom_point(
      position = position_jitter(width = 0.15),
      size = 3
    ) +
    
    scale_fill_manual(values = cfg$colors,
                      labels = c("luminal",
                                 "mucosal")) +
    
    scale_color_manual(values = cfg$colors,
                       labels = c("luminal",
                                  "mucosal")) +
    
    scale_x_discrete(labels = c("luminal",
                                "mucosal")) +
    
    theme_prism(base_size = 16) +
    
    labs(
      title = cfg$title,
      x = NULL,
      y = "Shannon index"
    ) +
    
    theme(
      axis.text.x = element_text(face = "bold"),
      legend.position = "none"
    ) +
    
    coord_cartesian(ylim = y_axis_limits) +
    
    geom_signif(
      comparisons = list(c("scraping","content")),
      annotations = p_label,
      y_position = cont_max + 0.01,
      tip_length = 0,
      size = 1.2,
      textsize = 8,
      color = "black"
    ) 
  ggsave(cfg$filename,
         p,
         width = 8,
         height = 8)
}



