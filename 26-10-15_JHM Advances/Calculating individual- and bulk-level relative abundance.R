library(tidyverse)
library(ggplot2)
library(RColorBrewer)
df <- read.csv(
  "User microbiome data.csv",
  check.names = FALSE
)
out_dir <- "User workstation pathway"
if (!dir.exists(out_dir)) dir.create(out_dir)
sample_cols <- grep("^(Cont|Direct|Indirect)\\d+$", colnames(df), value = TRUE)
levels_to_plot <- c("Phylum","Class","Order","Family","Genus","Species")
levels_to_plot <- levels_to_plot[levels_to_plot %in% colnames(df)]
make_palette <- function(tax_levels) {
  tax_names <- setdiff(tax_levels, c("ETC", "Unassigned"))
  n_colors  <- length(tax_names)
  base_cols_full <- c(
    "#1D70A8",
    "#FF7F00",
    "#4DAF1A",
    "#F0017D",
    "#6A3C9B",
    "#A65628",
    "#F781BF",
    "#000000",
    "#1B9E77",
    "#D95F02"
  )
  base_cols <- base_cols_full[seq_len(min(n_colors, length(base_cols_full)))]
  palette_vec <- setNames(base_cols, tax_names[seq_along(base_cols)])
  palette_colors <- c(
    palette_vec,
    "ETC"        = "#999999",
    "Unassigned" = "#CCCCCC"
  )
  palette_colors <- palette_colors[tax_levels]
  palette_colors
}
cliff_delta_safe <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }
  
  comparison_mat <- outer(x, y, "-")
  
  delta <- (
    sum(comparison_mat > 0) -
      sum(comparison_mat < 0)
  ) / (length(x) * length(y))
  
  as.numeric(delta)
}
for (lvl in levels_to_plot) {
  mdf <- df %>%
    select(all_of(lvl), all_of(sample_cols)) %>%
    pivot_longer(
      -all_of(lvl),
      names_to  = "Sample",
      values_to = "Count"
    ) %>%
    mutate(
      Count = as.numeric(Count),
      !!lvl := ifelse(
        is.na(.data[[lvl]]) | .data[[lvl]] == "" | .data[[lvl]] == "-",
        "Unassigned",
        .data[[lvl]]
      ),
      Group = case_when(
        grepl("^Cont",     Sample) ~ "Cont",
        grepl("^Direct",   Sample) ~ "Direct",
        grepl("^Indirect", Sample) ~ "Indirect",
        TRUE                        ~ NA_character_
      )
    ) %>%
    filter(!is.na(Count), !is.na(Group)) %>%
    group_by(Sample) %>%
    mutate(RelAbund = Count / sum(Count, na.rm = TRUE)) %>%
    ungroup()
  top_taxa <- mdf %>%
    group_by(.data[[lvl]]) %>%
    summarise(mean_abund = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    pull(.data[[lvl]]) %>%
    head(10)
  if (!"Unassigned" %in% top_taxa) {
    top_taxa <- c(top_taxa, "Unassigned")
  }
  mdf <- mdf %>%
    mutate(
      Tax_Filtered = ifelse(
        .data[[lvl]] %in% top_taxa,
        .data[[lvl]],
        "ETC"
      )
    )
  sample_df <- mdf %>%
    group_by(Sample, Group, Tax_Filtered) %>%
    summarise(RelAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(RelAbund = RelAbund / sum(RelAbund, na.rm = TRUE)) %>%
    ungroup()
  tax_order <- sample_df %>%
    group_by(Tax_Filtered) %>%
    summarise(mean_abund = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    pull(Tax_Filtered)
  sample_df$Tax_Filtered <- factor(sample_df$Tax_Filtered, levels = tax_order)
  sample_df$Sample <- factor(sample_df$Sample, levels = sort(unique(sample_df$Sample)))
  pal_colors <- make_palette(tax_order)
  p_sample <- ggplot(sample_df,
                     aes(x = Sample, y = RelAbund, fill = Tax_Filtered)) +
    geom_bar(
      stat     = "identity",
      position = "stack",
      color    = NA,
      width    = 0.9
    ) +
    scale_fill_manual(values = pal_colors) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.text.x   = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y   = element_text(size = 10),
      axis.title.x  = element_text(size = 24, face = "bold"),
      axis.title.y  = element_text(size = 24, face = "bold"),
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title  = element_text(size = 20, face = "bold"),
      legend.text   = element_text(size = 10)
    ) +
    labs(
      title = paste(lvl, "relative abundance by sample"),
      x     = "Sample",
      y     = "Relative abundance",
      fill  = lvl
    )
  ggsave(
    filename = file.path(out_dir, paste0("Relative_Abundance_Sample_", lvl, ".png")),
    plot     = p_sample,
    width    = 10,
    height   = 5,
    dpi      = 600
  )
  write.csv(
    sample_df,
    file      = file.path(out_dir, paste0("Relative_Abundance_Sample_", lvl, ".csv")),
    row.names = FALSE
  )
  group_df <- sample_df %>%
    group_by(Group, Tax_Filtered) %>%
    summarise(mean_rel = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(
      rel_in_group = mean_rel / sum(mean_rel, na.rm = TRUE)
    ) %>%
    ungroup()
  group_df$Tax_Filtered <- factor(group_df$Tax_Filtered, levels = tax_order)
  group_df$Group <- factor(group_df$Group, levels = c("Cont", "Indirect", "Direct"))
  p_group <- ggplot(group_df,
                    aes(x = Group, y = rel_in_group, fill = Tax_Filtered)) +
    geom_bar(
      stat     = "identity",
      position = "stack",
      color    = NA,
      width    = 0.9
    ) +
    scale_fill_manual(values = pal_colors) +
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.text.x   = element_text(size = 18, angle = 45, hjust = 1),
      axis.text.y   = element_text(size = 18),
      axis.title.x  = element_text(size = 24, face = "bold"),
      axis.title.y  = element_text(size = 24, face = "bold"),
      plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title  = element_text(size = 25, face = "bold"),
      legend.text   = element_text(size = 15)
    ) +
    labs(
      title = paste(lvl, "relative abundance within group (sum = 1)"),
      x     = "Group",
      y     = "Relative abundance within group (sum = 1)",
      fill  = lvl
    )
  ggsave(
    filename = file.path(out_dir, paste0("Relative_Abundance_GroupNorm_", lvl, ".png")),
    plot     = p_group,
    width    = 8,
    height   = 12,
    dpi      = 600
  )
  write.csv(
    group_df,
    file      = file.path(out_dir, paste0("Relative_Abundance_GroupNorm_", lvl, ".csv")),
    row.names = FALSE
  )
}