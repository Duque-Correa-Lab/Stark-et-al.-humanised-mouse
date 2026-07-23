# Shannon Index in uncultured metascrapes, and in metascarapes after aerobic and anaerobic culture
#
# ============================
# 1. Load libraries
# ============================
library(tidyverse)
library(vegan)
library(ggprism)      
library(ggsignif)     
library(dunn.test)    
library(dplyr)
library(ggplot2)

setwd("/Your/path/here") #add correct path

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

# ============================
# 3. Keep species only
# ============================
otu_scraping_sp <- otu_scraping %>%
  filter(grepl("s__", clade_name), !grepl("t__", clade_name)) %>%
  distinct(clade_name, .keep_all = TRUE)

otu_mat_scraping <- otu_scraping_sp %>%
  column_to_rownames("clade_name") %>%
  as.matrix()

# transpose
otu_mat_scraping <- t(otu_mat_scraping)

# clean sample names
rownames(otu_mat_scraping) <- sub("_results_gtdb", "", rownames(otu_mat_scraping))


# ============================
# 4. Match sample order between read-count file and OTU table
# ============================

# Define shared samples
shared_samples <- intersect(
  rownames(otu_mat_scraping),
  seq_scraping$`Run-Accession`
)

seq_scraping_filt <- seq_scraping %>%
  filter(`Run-Accession` %in% shared_samples) %>%
  arrange(match(`Run-Accession`, shared_samples))

# Final order alignment
otu_mat_scraping <- otu_mat_scraping[seq_scraping_filt$`Run-Accession`, ]

# Verify alignment safely
stopifnot(identical(rownames(otu_mat_scraping), seq_scraping_filt$`Run-Accession`))

# ============================
# 5. Convert relative abundance → counts
# ============================
otu_mat_prop <- otu_mat_scraping / 100

otu_mat_counts <- otu_mat_prop
for(i in 1:nrow(otu_mat_counts)){
  otu_mat_counts[i, ] <- round(otu_mat_prop[i, ] * seq_scraping_filt$`Read-count`[i])
}

otu_mat_counts[is.na(otu_mat_counts)] <- 0
mode(otu_mat_counts) <- "integer"


# ============================
# 6. Rarefy separately per annot1 group
# ============================
meta_ordered <- meta_scraping %>%
  filter(run_accession %in% rownames(otu_mat_counts)) %>%
  arrange(match(run_accession, rownames(otu_mat_counts)))

otu_mat_counts <- otu_mat_counts[meta_ordered$run_accession, ]

# Verify alignment (same check used in step 4)
stopifnot(identical(rownames(otu_mat_counts), meta_ordered$run_accession))

set.seed(123)

rarefied_list <- list()

for(group in unique(meta_ordered$annot1)) {
  
  samples_group <- meta_ordered$run_accession[meta_ordered$annot1 == group]
  
  mat_group <- otu_mat_counts[samples_group, , drop = FALSE]
  
  mat_group <- mat_group[rowSums(mat_group) > 0, , drop = FALSE]
  
  # Skip groups with too few samples to rarefy meaningfully (same guard
  # used in the MS-vs-LC script) - avoids min()/rrarefy() erroring out on
  # an empty or single-row matrix
  if(nrow(mat_group) < 2) next
  
  min_depth <- min(rowSums(mat_group))
  
  rarefied_list[[group]] <- vegan::rrarefy(mat_group, sample = min_depth)
}

otu_rarefied_scraping <- do.call(rbind, rarefied_list)



# ============================
# 7. Shannon diversity
# ============================
shannon_scraping <- diversity(otu_rarefied_scraping, index = "shannon")

shannon_df_scraping <- tibble(
  run_accession = rownames(otu_rarefied_scraping),
  Shannon = shannon_scraping
)

# ============================
# 8. Merge metadata 
# ============================
shannon_meta_scraping <- shannon_df_scraping %>%
  left_join(meta_scraping, by = "run_accession")

# ============================
# 9. Plot: alpha diversity within each annot1 group comparing annot2 
# ============================

plot_configs <- list(
  
  list(
    filter_vals = c("WT MS Uncultured", "WT MS Aerobic", "WT MS Anaerobic"),
    colors = c("WT MS Uncultured" = "gray12",
               "WT MS Aerobic" = "gray12",
               "WT MS Anaerobic" = "gray12"),
    shapes = c(15, 16, 17),
    filename = "Shannon_depthAdj_WT_MS_all_conditions.pdf",
    title = "WT MS: Culture conditions"
  ),
  
  list(
    filter_vals = c("D2 HMA MS Uncultured", "D2 HMA MS Aerobic", "D2 HMA MS Anaerobic"),
    colors = c("D2 HMA MS Uncultured" = "steelblue4",
               "D2 HMA MS Aerobic" = "steelblue4",
               "D2 HMA MS Anaerobic" = "steelblue4"),
    shapes = c(15, 16, 17),
    filename = "Shannon_depthAdj_D2HMA_MS_all_conditions.pdf",
    title = "D2 HMA MS: Culture conditions "
  ),
  
  list(
    filter_vals = c("D7 HMA MS Uncultured", "D7 HMA MS Aerobic", "D7 HMA MS Anaerobic"),
    colors = c("D7 HMA MS Uncultured" = "indianred4",
               "D7 HMA MS Aerobic" = "indianred4",
               "D7 HMA MS Anaerobic" = "indianred4"),
    shapes = c(15, 16, 17),
    filename = "Shannon_depthAdj_D7HMA_MS_all_conditions.pdf",
    title = "D7 HMA MS: Culture conditions"
  )
  
)

# ============================
# Shared y-axis range across all panels (so plots are directly comparable
# side by side for publication)
# ============================
overall_ymin <- min(shannon_meta_scraping$Shannon, na.rm = TRUE)
overall_ymax <- max(shannon_meta_scraping$Shannon, na.rm = TRUE)
overall_yrange <- overall_ymax - overall_ymin

# Headroom reserved for the worst case: 3 groups per panel -> up to 3
# stacked significance brackets. Fixed regardless of how many are actually
# significant in any given panel, so the axis is identical across panels.
max_possible_comparisons <- 3
y_axis_limits <- c(
  overall_ymin - overall_yrange * 0.05,
  overall_ymax + overall_yrange * (0.08 + 0.15 * max_possible_comparisons + 0.05)
)

for(cfg in plot_configs) {
  
  df <- shannon_meta_scraping %>%
    filter(annot %in% cfg$filter_vals) %>%
    mutate(annot = factor(annot, levels = cfg$filter_vals),
           shape_group = factor(annot2,
                                levels = c("Uncultured", "Aerobic", "Anaerobic")))
  
  # ============================
  # Kruskal-Wallis
  # ============================
  
  kw <- kruskal.test(Shannon ~ annot, data = df)
  
  print(cfg$title)
  print(kw)
  
  # ============================
  # Dunn test
  # ============================
  
  dunn <- dunn.test(df$Shannon, df$annot, kw = TRUE, label = TRUE)
  print(dunn)
  
  # Extract significant comparisons
  sig_comparisons <- list()
  sig_labels <- c()
  
  for(i in seq_along(dunn$comparisons)) {
    
    if(dunn$P.adjusted[i] < 0.05) {
      
      groups <- strsplit(dunn$comparisons[i], " - ")[[1]]
      
      sig_comparisons[[length(sig_comparisons)+1]] <- groups
      
      sig_labels <- c(sig_labels, p_to_stars(dunn$P.adjusted[i]))
    }
  }
  
  # ============================
  # Plot
  # ============================
  
  p <- ggplot(df,
              aes(x = annot,
                  y = Shannon,
                  fill = annot,
                  color = annot)) +
    
    geom_boxplot(
      color = "black",
      width = 0.55,
      coef = 0,
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
      aes(shape = shape_group),
      position = position_jitter(width = 0.15),
      size = 3
    ) +
    
    scale_fill_manual(values = cfg$colors) +
    scale_color_manual(values = cfg$colors) +
    scale_shape_manual(values = cfg$shapes) +
    
    # Rename x-axis to simple names
    scale_x_discrete(labels = setNames(
      c("uncultured", "aerobic", "anaerobic"),
      cfg$filter_vals
    )) +
    
    theme_prism(base_size = 16) +
    
    labs(
      title = cfg$title,
      x = NULL,
      y = "Shannon-index"
    ) +
    
    theme(
      axis.text.x = element_text(angle = 0,
                                 hjust = 0.5,
                                 face = "bold"),
      legend.position = "none"
    ) +
    
    coord_cartesian(ylim = y_axis_limits)
  
  # ============================
  # Add significance bars
  # ============================
  
  if(length(sig_comparisons) > 0) {
    
    ymax <- max(df$Shannon)
    n_sig <- length(sig_comparisons)
    
    # Bracket spacing based on the shared overall range (not this panel's
    # own range) so bracket heights line up consistently across panels
    bracket_step <- overall_yrange * 0.15
    y_positions <- ymax + overall_yrange*0.08 + bracket_step * (0:(n_sig - 1))
    
    p <- p + geom_signif(
      comparisons = sig_comparisons,
      annotations = sig_labels,
      y_position = y_positions,
      tip_length = 0.01,
      size = 1.2,
      color = "black",
      textsize = 5
    )
  }
  
  # ============================
  # Save
  # ============================
  
  ggsave(cfg$filename,
         p,
         width = 6,
         height = 6)
  
}

