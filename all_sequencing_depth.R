library(dplyr)
library(readr)
library(ggplot2)
library(ggprism)
library(rstatix)
library(ggpubr)

#Set dir
setwd("/Users/your/path/here")
meta_dir <- "/Users/your/path/here"

# Read each file
seq_content  <- read_csv(file.path(meta_dir, "read_counts_content.csv"))
seq_scraping <- read_csv(file.path(meta_dir, "read_counts_scraping.csv"))
seq_wc       <- read_csv(file.path(meta_dir, "read_counts_whole_caecum.csv"))

# Combine all datasets
seq_all <- bind_rows(seq_content, seq_scraping, seq_wc)

#Read in metadata
meta_content <- read_tsv(file.path(meta_dir, "AerobevsAnaerobe_content_metadata.txt")) %>%
  dplyr::select(run_accession, annot)

meta_scraping <- read_tsv(file.path(meta_dir, "AerobevsAnaerobe_scraping_metadata.txt")) %>%
  dplyr::select(run_accession, annot)

meta_wc <- read_tsv(file.path(meta_dir, "metadata_WC_HMAvsWT.txt")) %>%
  dplyr::select(run_accession, annot)


# Check for duplicated run_accession
dup_content  <- meta_content  %>% count(run_accession) %>% filter(n > 1)
dup_scraping <- meta_scraping %>% count(run_accession) %>% filter(n > 1)
dup_wc       <- meta_wc       %>% count(run_accession) %>% filter(n > 1)

# dup_content
dup_content
dup_scraping
dup_wc


# Stack all metadata tables into one long table
meta_merged <- bind_rows(meta_content, meta_scraping, meta_wc)


#Merge metadata and seq_annot by run_accession
seq_annot <- seq_all %>%
  left_join(meta_merged, by = c('Run-Accession' = "run_accession"))

#Flag samples not in metadata
seq_annot <- seq_annot %>%
  mutate(
    in_metadata = !is.na(annot)
  )

# Check missing
seq_annot %>%
  filter(!in_metadata)

# Create ordered factor columns
seq_annot <- seq_annot %>%
  mutate(
    # Host status order
    host = case_when(
      grepl("^WT", annot) ~ "WT",
      grepl("^D2 HMA", annot) ~ "D2 HMA",
      grepl("^D7 HMA", annot) ~ "D7 HMA",
      TRUE ~ NA_character_
    ),
    
    host = factor(host, levels = c("WT", "D2 HMA", "D7 HMA")),
    
    # Sample type order
    material = case_when(
      grepl("MS", annot) ~ "MS",
      grepl("LC", annot) ~ "LC",
      TRUE ~ NA_character_
    ),
    
    material = factor(material, levels = c("MS", "LC")),
    
    # Culture condition order
    culture = case_when(
      grepl("Uncultured", annot) ~ "Uncultured",
      grepl("Aerobic", annot) ~ "Aerobic",
      grepl("Anaerobic", annot) ~ "Anaerobic",
      TRUE ~ NA_character_
    ),
    
    culture = factor(culture, levels = c("Uncultured", "Aerobic", "Anaerobic"))
  )


# Ensure seq_annot has a combined annot column for plotting
seq_annot <- seq_annot %>%
  mutate(annot = paste(host, material, culture, sep = " "))

# Function to reorder annot by culture for MS-only plots
reorder_by_culture <- function(df) {
  df %>%
    mutate(
      annot = factor(
        annot,
        levels = paste(unique(df$host), unique(df$material), c("Uncultured", "Aerobic", "Anaerobic"))
      )
    )
}



# ==============================================================================
# Adaptive Stat Function (creates significance star labels)
# ==============================================================================
add_stats <- function(df) {
  df$annot <- factor(df$annot)
  unique_groups <- length(levels(df$annot))
  
  if (unique_groups < 2) {
    message("Fewer than 2 active groups found. Skipping statistical testing.")
    return(NULL)
  }
  
  if (unique_groups == 2) {
    message("Processing 2 groups. Executing Wilcoxon Test...")
    stat.test <- df %>%
      wilcox_test(`Read-count` ~ annot) %>%
      add_significance(p.col = "p", output.col = "p.adj.signif") %>%
      mutate(p.adj = p) %>%
      add_xy_position(x = "annot")
  } else {
    message(paste("Processing", unique_groups, "groups. Executing Kruskal-Wallis + Dunn..."))
    kw <- df %>% kruskal_test(`Read-count` ~ annot)
    print(kw)
    stat.test <- df %>%
      dunn_test(`Read-count` ~ annot, p.adjust.method = "bonferroni") %>%
      add_xy_position(x = "annot")
  }
  print(stat.test, width = Inf)
  stat.test
}

# Draws brackets manually in log10 space -- bypasses an earlier broken stat_pvalue_manual/log-scale interaction
add_bracket_layers <- function(p, stat.test, df, y_col = "Read-count") {
  if (is.null(stat.test)) return(p)
  sig <- stat.test %>% filter(p.adj.signif != "ns")
  if (nrow(sig) == 0) return(p)
  
  y_max_log <- log10(max(df[[y_col]], na.rm = TRUE))
  sig <- sig %>%
    mutate(y.log = y_max_log + 0.06 + 0.08 * (row_number() - 1), tip = 0.015)
  
  for (i in seq_len(nrow(sig))) {
    r  <- sig[i, ]
    y  <- 10^r$y.log
    y0 <- 10^(r$y.log - r$tip)
    p <- p +
      annotate("segment", x = r$xmin, xend = r$xmin, y = y0, yend = y) +
      annotate("segment", x = r$xmax, xend = r$xmax, y = y0, yend = y) +
      annotate("segment", x = r$xmin, xend = r$xmax, y = y,  yend = y) +
      annotate("text", x = (r$xmin + r$xmax) / 2, y = 10^(r$y.log + 0.03),
               label = r$p.adj.signif, size = 5)
  }
  p
}


# ==============================================================================
# Reusable Plotting Function
# ==============================================================================
plot_reads = function(df, title,
                      colors,
                      shapes = NULL,
                      alpha_box = 0.5,
                      median_width = 0.3,
                      y_limits = NULL) {
  
  # Compute stats adaptively (Returns Wilcoxon structure or Dunn structure)
  stat.test = add_stats(df)
  
  p = ggplot(df, aes(x = annot, y = `Read-count`)) +
    geom_boxplot(
      aes(fill = annot),
      color = "black",
      outlier.shape = NA,
      alpha = alpha_box,
      width = 0.6,
      linewidth = 0.8,
      coef = 0
    ) +
    stat_summary(
      aes(color = annot),
      fun = median,
      geom = "crossbar",
      width = 0.6,
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    geom_jitter(
      aes(color = annot, shape = annot),
      width = 0.2,
      size = 4,
      alpha = 0.9
    ) +
    scale_y_log10(limits = c(1e5, 2e7)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(
      x = NULL,
      y = "Number of reads",
      title = title
    ) +
    theme_prism(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none"
    )
  
  # Plot significance bar overlays if tests ran successfully and found results
  if (!is.null(stat.test) && nrow(stat.test) > 0) {
    p <- add_bracket_layers(p, stat.test, df)
  }
  
  if (!is.null(shapes)) {
    p = p + scale_shape_manual(values = shapes)
  } else {
    p = p + scale_shape_manual(values = rep(16, length(unique(df$annot))))
  }
  
  return(p)
}

# ==============================================================================
# Generate plots for comparisons
# ==============================================================================

# 1. WT vs D2 HMA vs D7 HMA (3 Groups -> Runs Kruskal-Wallis Automatically)
df_1 = seq_annot %>%
  filter(annot %in% c("WT NA NA", "D2 HMA NA NA", "D7 HMA NA NA")) %>%
  mutate(annot = factor(annot, levels = c("WT NA NA", "D2 HMA NA NA", "D7 HMA NA NA")))

colors_1 = c(
  "WT NA NA"     = "gray45",
  "D2 HMA NA NA" = "steelblue",
  "D7 HMA NA NA" = "indianred"
)

p1 = plot_reads(df_1, title = "Whole caecum", colors = colors_1, shapes = c(16, 16, 16)) +
  scale_x_discrete(labels = c("WT NA NA" = "WT", "D2 HMA NA NA" = "D2 HMA", "D7 HMA NA NA" = "D7 HMA"))
p1

# 2. WT only: MS vs LC Uncultured (2 Groups -> Runs Wilcoxon Automatically)
df_2 = seq_annot %>% filter(annot %in% c("WT MS Uncultured", "WT LC Uncultured"))
p2 = plot_reads(df_2, title = "WT", colors = c("WT MS Uncultured" = "gray12", "WT LC Uncultured" = "gray40"), shapes = c(16, 16)) +
  scale_x_discrete(labels = c("WT MS Uncultured" = "mucosal", "WT LC Uncultured" = "luminal"))
p2

# 3. D2 HMA only: MS vs LC Uncultured (2 Groups -> Runs Wilcoxon Automatically)
df_3 = seq_annot %>% filter(annot %in% c("D2 HMA MS Uncultured", "D2 HMA LC Uncultured"))
p3 = plot_reads(df_3, title = "D2 HMA", colors = c("D2 HMA MS Uncultured" = "steelblue4", "D2 HMA LC Uncultured" = "steelblue2"), shapes = c(16, 16)) +
  scale_x_discrete(labels = c("D2 HMA MS Uncultured" = "mucosal", "D2 HMA LC Uncultured" = "luminal"))
p3

# 4. D7 HMA only: MS vs LC Uncultured (2 Groups -> Runs Wilcoxon Automatically)
df_4 = seq_annot %>% filter(annot %in% c("D7 HMA MS Uncultured", "D7 HMA LC Uncultured"))
p4 = plot_reads(df_4, title = "D7 HMA", colors = c("D7 HMA MS Uncultured" = "#9E2A2F", "D7 HMA LC Uncultured" = "#FA8"), shapes = c(16, 16)) +
  scale_x_discrete(labels = c("D7 HMA MS Uncultured" = "mucosal", "D7 HMA LC Uncultured" = "luminal"))
p4

# 5. WT MS: Uncultured → Aerobic → Anaerobic (3 Groups -> Runs Kruskal-Wallis Automatically)
df_5 = seq_annot %>% filter(host == "WT", material == "MS") %>% reorder_by_culture()
p5 = plot_reads(df_5, title = "WT mucosal scraping", colors = c("WT MS Uncultured" = "gray12", "WT MS Aerobic" = "gray12", "WT MS Anaerobic" = "gray12"), shapes = c(15, 16, 17)) +
  scale_x_discrete(labels = c("WT MS Uncultured" = "uncultured", "WT MS Aerobic" = "aerobic", "WT MS Anaerobic" = "anaerobic"))
p5

# 6. D2 HMA MS: Uncultured → Aerobic → Anaerobic (3 Groups -> Runs Kruskal-Wallis Automatically)
df_6 = seq_annot %>% filter(host == "D2 HMA", material == "MS") %>% reorder_by_culture()
p6 = plot_reads(df_6, title = "D2 HMA mucosal scrapings", colors = c("D2 HMA MS Uncultured" = "steelblue4", "D2 HMA MS Aerobic" = "steelblue4", "D2 HMA MS Anaerobic" = "steelblue4"), shapes = c(15, 16, 17)) +
  scale_x_discrete(labels = c("D2 HMA MS Uncultured" = "uncultured", "D2 HMA MS Aerobic" = "aerobic", "D2 HMA MS Anaerobic" = "anaerobic"))
p6

# 7. D7 HMA MS: Uncultured → Aerobic → Anaerobic (3 Groups -> Runs Kruskal-Wallis Automatically)
df_7 = seq_annot %>% filter(host == "D7 HMA", material == "MS") %>% reorder_by_culture()
p7 = plot_reads(df_7, title = "D7 HMA mucosal scrapings", colors = c("D7 HMA MS Uncultured" = "#9E2A2F", "D7 HMA MS Aerobic" = "#9E2A2F", "D7 HMA MS Anaerobic" = "#9E2A2F"), shapes = c(15, 16, 17)) +
  scale_x_discrete(labels = c("D7 HMA MS Uncultured" = "uncultured", "D7 HMA MS Aerobic" = "aerobic", "D7 HMA MS Anaerobic" = "anaerobic"))
p7

# Create a list of all plots with descriptive names
plots_list <- list(
  "sequencing_depth_WT_vs_D2_vs_D7"   = p1,
  "sequencing_depth_MSvsLC_WT"        = p2,
  "sequencing_depth_MSvsLC_D2HMA"     = p3,
  "sequencing_depth_MSvsLC_D7HMA"     = p4,
  "sequencing_depth_AervsAnaer_MS_WT" = p5,
  "sequencing_depth_AervsAnaer_MS_D2HMA" = p6,
  "sequencing_depth_AervsAnaer_MS_D7HMA" = p7
)

# Set output directory for saved plots
output_dir <- file.path(meta_dir, "Sequencing_depth_plots_updated")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Loop over plots and save each
for (name in names(plots_list)) {
  ggsave(
    filename = paste0(name, ".pdf"),
    plot = plots_list[[name]],
    path = output_dir,
    width = 8,
    height = 6,
    dpi = 300
  )
}


