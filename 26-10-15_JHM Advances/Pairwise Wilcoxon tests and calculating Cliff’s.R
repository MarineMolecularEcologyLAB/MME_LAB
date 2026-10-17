library(tidyverse)
library(ggplot2)
library(effsize)

df <- read.csv(
  "User microbiome data.csv ",
  check.names = FALSE
)

out_dir <- "User workstation pathway"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

sample_cols <- grep(
  "^(Cont|Direct|Indirect)\\d+$",
  colnames(df),
  value = TRUE
)

levels_to_plot <- c(
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

levels_to_plot <- levels_to_plot[
  levels_to_plot %in% colnames(df)
]

reference_group <- "Cont"
compare_groups <- c("Direct", "Indirect")

do_sample_plot <- TRUE
do_mean_plot <- TRUE

for (lvl in levels_to_plot) {
  
  mdf_lvl <- df %>%
    dplyr::select(
      all_of(lvl),
      all_of(sample_cols)
    ) %>%
    pivot_longer(
      cols = -all_of(lvl),
      names_to = "Sample",
      values_to = "Count"
    ) %>%
    mutate(
      Count = suppressWarnings(as.numeric(Count)),
      !!lvl := ifelse(
        is.na(.data[[lvl]]) |
          .data[[lvl]] == "" |
          .data[[lvl]] == "-",
        "Unassigned",
        as.character(.data[[lvl]])
      ),
      Group = case_when(
        grepl("^Cont", Sample) ~ "Cont",
        grepl("^Direct", Sample) ~ "Direct",
        grepl("^Indirect", Sample) ~ "Indirect",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(Count),
      !is.na(Group)
    ) %>%
    group_by(Sample) %>%
    mutate(
      sample_total = sum(Count, na.rm = TRUE),
      RelAbund = ifelse(
        sample_total > 0,
        Count / sample_total,
        0
      )
    ) %>%
    ungroup() %>%
    dplyr::select(-sample_total)
  
  res_delta <- expand.grid(
    Taxon = unique(mdf_lvl[[lvl]]),
    Group2 = compare_groups,
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      ref_vals = list(
        mdf_lvl$RelAbund[
          mdf_lvl[[lvl]] == Taxon &
            mdf_lvl$Group == reference_group
        ]
      ),
      cmp_vals = list(
        mdf_lvl$RelAbund[
          mdf_lvl[[lvl]] == Taxon &
            mdf_lvl$Group == Group2
        ]
      ),
      n_ref = length(unlist(ref_vals)),
      n_cmp = length(unlist(cmp_vals)),
      p = if (n_ref >= 2 && n_cmp >= 2) {
        tryCatch(
          suppressWarnings(
            wilcox.test(
              unlist(ref_vals),
              unlist(cmp_vals),
              exact = FALSE
            )$p.value
          ),
          error = function(e) NA_real_
        )
      } else {
        NA_real_
      },
      delta = if (n_ref >= 2 && n_cmp >= 2) {
        tryCatch(
          as.numeric(
            suppressWarnings(
              effsize::cliff.delta(
                unlist(cmp_vals),
                unlist(ref_vals)
              )$estimate
            )
          ),
          error = function(e) NA_real_
        )
      } else {
        NA_real_
      },
      mean_abund = mean(
        c(
          unlist(ref_vals),
          unlist(cmp_vals)
        ),
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    mutate(
      sig = !is.na(p) & p < 0.1,
      dir = case_when(
        is.na(delta) ~ NA_character_,
        delta >= 0 ~ "Positive association",
        delta < 0 ~ "Negative association"
      ),
      size_v = log10(
        pmax(mean_abund, 1e-6)
      )
    )
  
  write.csv(
    res_delta %>%
      dplyr::select(-ref_vals, -cmp_vals),
    file.path(
      out_dir,
      paste0(
        "Triangle_Wilcoxon_CliffDelta_",
        lvl,
        "_vs_",
        reference_group,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  sig_pairs <- res_delta %>%
    filter(
      sig,
      is.finite(delta),
      is.finite(size_v),
      !is.na(dir)
    ) %>%
    dplyr::select(
      Taxon,
      Group2,
      delta,
      dir,
      size_v
    )
  
  if (nrow(sig_pairs) == 0) {
    next
  }
  
  plot_df_sample <- sig_pairs %>%
    inner_join(
      mdf_lvl %>%
        filter(Group %in% compare_groups) %>%
        dplyr::select(
          Taxon = all_of(lvl),
          Sample,
          Group,
          RelAbund
        ),
      by = c(
        "Taxon",
        "Group2" = "Group"
      )
    ) %>%
    mutate(
      Sample = factor(
        Sample,
        levels = sort(unique(Sample))
      ),
      size_cat = cut(
        RelAbund,
        breaks = c(
          -Inf,
          1e-4,
          1e-3,
          1e-2,
          1e-1,
          Inf
        ),
        labels = c(
          "very low",
          "low",
          "mid",
          "high",
          "very high"
        ),
        right = FALSE
      )
    )
  
  plot_df_mean <- plot_df_sample %>%
    group_by(
      Taxon,
      Group2
    ) %>%
    summarise(
      RelAbund = mean(
        RelAbund,
        na.rm = TRUE
      ),
      delta = delta[1],
      dir = dir[1],
      .groups = "drop"
    ) %>%
    mutate(
      Group2 = factor(
        Group2,
        levels = compare_groups
      )
    )
  
  if (do_sample_plot && nrow(plot_df_sample) > 0) {
    
    p_tri_sample <- ggplot(
      plot_df_sample,
      aes(
        x = Taxon,
        y = Sample
      )
    ) +
      geom_point(
        aes(
          shape = dir,
          fill = delta,
          size = size_cat
        ),
        color = "black",
        stroke = 0.3,
        alpha = 0.8
      ) +
      scale_shape_manual(
        values = c(
          "Positive association" = 24,
          "Negative association" = 25
        ),
        name = "Direction"
      ) +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        limits = c(-1, 1),
        oob = scales::squish,
        name = "Effect size"
      ) +
      scale_size_manual(
        name = "Relative abundance",
        values = c(
          "very low" = 6,
          "low" = 10,
          "mid" = 14,
          "high" = 18,
          "very high" = 22
        ),
        drop = FALSE,
        guide = guide_legend(
          override.aes = list(
            shape = 24,
            fill = "black",
            colour = "black"
          )
        )
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(
          size = 14,
          face = "bold"
        ),
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        legend.title = element_text(
          size = 12,
          face = "bold"
        )
      ) +
      labs(
        title = paste(
          lvl,
          "significant changes vs",
          reference_group,
          "(sample-level)"
        ),
        x = lvl,
        y = "Sample"
      )
    
    sample_png <- file.path(
      out_dir,
      paste0(
        "Triangle_SAMPLE_",
        lvl,
        "_vs_",
        reference_group,
        ".png"
      )
    )
    
    tryCatch(
      ggsave(
        filename = sample_png,
        plot = p_tri_sample,
        device = "png",
        width = 14,
        height = 10,
        units = "in",
        dpi = 600,
        limitsize = FALSE
      ),
      error = function(e) {
        png(
          filename = sample_png,
          width = 8400,
          height = 6000,
          res = 600
        )
        print(p_tri_sample)
        dev.off()
      }
    )
    
    write.csv(
      plot_df_sample,
      file.path(
        out_dir,
        paste0(
          "Triangle_SAMPLE_TABLE_",
          lvl,
          "_vs_",
          reference_group,
          ".csv"
        )
      ),
      row.names = FALSE
    )
  }
  
  if (do_mean_plot && nrow(plot_df_mean) > 0) {
    
    p_tri_mean <- ggplot(
      plot_df_mean,
      aes(
        x = Taxon,
        y = Group2
      )
    ) +
      geom_point(
        aes(
          shape = dir,
          fill = delta,
          size = RelAbund
        ),
        color = "black",
        stroke = 0.7
      ) +
      scale_shape_manual(
        values = c(
          "Positive association" = 24,
          "Negative association" = 25
        ),
        name = "Direction"
      ) +
      scale_fill_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        limits = c(-1, 1),
        oob = scales::squish,
        name = "Effect size (Cliff's delta)"
      ) +
      scale_size(
        range = c(3, 7),
        name = "Mean relative abundance"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(
          size = 14,
          face = "bold"
        ),
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        legend.title = element_text(
          size = 12,
          face = "bold"
        )
      ) +
      labs(
        title = paste(
          lvl,
          "significant changes vs",
          reference_group,
          "(group mean)"
        ),
        x = lvl,
        y = "Group"
      )
    
    mean_png <- file.path(
      out_dir,
      paste0(
        "Triangle_MEAN_",
        lvl,
        "_vs_",
        reference_group,
        ".png"
      )
    )
    
    tryCatch(
      ggsave(
        filename = mean_png,
        plot = p_tri_mean,
        device = "png",
        width = 14,
        height = 10,
        units = "in",
        dpi = 600,
        limitsize = FALSE
      ),
      error = function(e) {
        png(
          filename = mean_png,
          width = 8400,
          height = 6000,
          res = 600
        )
        print(p_tri_mean)
        dev.off()
      }
    )
    
    write.csv(
      plot_df_mean,
      file.path(
        out_dir,
        paste0(
          "Triangle_MEAN_TABLE_",
          lvl,
          "_vs_",
          reference_group,
          ".csv"
        )
      ),
      row.names = FALSE
    )
  }
}