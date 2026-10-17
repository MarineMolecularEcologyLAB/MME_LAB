library(tidyverse)
library(ggplot2)
library(RColorBrewer)

## 기본 세팅 ----------------------------------------------------------
df <- read.csv(
  "C:/Users/Marine/Desktop/HM_Microbiome/Microbiome_Merged_All_withTaxonomy_251123.csv",
  check.names = FALSE
)

out_dir <- "C:/Users/Marine/Desktop/251123_Relative A"
if (!dir.exists(out_dir)) dir.create(out_dir)

sample_cols <- grep("^(Cont|Direct|Indirect)\\d+$", colnames(df), value = TRUE)  # 수정됨

levels_to_plot <- c("Phylum","Class","Order","Family","Genus","Species")
levels_to_plot <- levels_to_plot[levels_to_plot %in% colnames(df)]

print(paste("샘플 열 개수:", length(sample_cols)))
print("샘플 열 이름들:")
print(sample_cols)


## 공통 함수: 색 팔레트 ------------------------------------------------

make_palette <- function(tax_levels) {
  tax_names <- setdiff(tax_levels, c("ETC", "Unassigned"))
  n_colors  <- length(tax_names)
  
  base_cols_full <- c(
    "#1D70A8",  # blue-ish
    "#FF7F00",  # orange
    "#4DAF1A",  # green
    "#F0017D",  # magenta
    "#6A3C9B",  # purple
    "#A65628",  # brown
    "#F781BF",  # pink
    "#000000",  # black
    "#1B9E77",  # teal
    "#D95F02"   # dark orange
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

## 메인 루프 ----------------------------------------------------------

for (lvl in levels_to_plot) {
  cat("\n", lvl, "수준 결과...\n")
  
  ## 1) long-format + Sample별 RelAbund ------------------------------
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
  
  ## 2) 상위 taxa + ETC -----------------------------------------------
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
  
  ## 3) 샘플 barplot용 데이터 -----------------------------------------
  sample_df <- mdf %>%
    group_by(Sample, Group, Tax_Filtered) %>%
    summarise(RelAbund = sum(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    group_by(Sample) %>%
    mutate(RelAbund = RelAbund / sum(RelAbund, na.rm = TRUE)) %>%  # 샘플 내 합 1
    ungroup()
  
  tax_order <- sample_df %>%
    group_by(Tax_Filtered) %>%
    summarise(mean_abund = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abund)) %>%
    pull(Tax_Filtered)
  
  sample_df$Tax_Filtered <- factor(sample_df$Tax_Filtered, levels = tax_order)
  sample_df$Sample <- factor(sample_df$Sample, levels = sort(unique(sample_df$Sample)))
  
  pal_colors <- make_palette(tax_order)
  
  ## 3-1) 모든 샘플 (n수)
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
      panel.grid.major = element_blank(),   # 주요 격자선 제거
      panel.grid.minor = element_blank(),   # 보조 격자선 제거
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
  
  ## 4) 그룹 barplot용 데이터 (그룹 내 합 = 1) ------------------------
  group_df <- sample_df %>%
    group_by(Group, Tax_Filtered) %>%
    summarise(mean_rel = mean(RelAbund, na.rm = TRUE), .groups = "drop") %>%
    group_by(Group) %>%
    mutate(
      rel_in_group = mean_rel / sum(mean_rel, na.rm = TRUE)  # 각 그룹 막대 합 1
    ) %>%
    ungroup()
  
  group_df$Tax_Filtered <- factor(group_df$Tax_Filtered, levels = tax_order)
  group_df$Group <- factor(group_df$Group, levels = c("Cont", "Indirect", "Direct"))
  
  ## 4-1) 그룹 barplot -------------------------------------------------
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
  


  #########################
  ########################
  #######################

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(effsize)   # cliff.delta

## ===== 기본 세팅 =====
df <- read.csv("C:/Users/Marine/Desktop/HM_Microbiome/Microbiome_Merged_All_withTaxonomy_251123.csv",
               check.names = FALSE)

out_dir <- "C:/Users/Marine/Desktop/251123_Relative A"
if (!dir.exists(out_dir)) dir.create(out_dir)

sample_cols <- grep("^(Cont|Direct|Indirect)\\d+$", colnames(df), value = TRUE)
levels_to_plot <- c("Phylum","Class","Order","Family","Genus","Species")
levels_to_plot <- levels_to_plot[levels_to_plot %in% colnames(df)]

reference_group <- "Cont"
compare_groups  <- c("Direct", "Indirect")
do_sample_plot <- TRUE
do_mean_plot   <- TRUE

## ===== level 루프 시작 =====
for (lvl in levels_to_plot) {
  
  cat("\n[TRIANGLE] ", lvl, " 수준 Wilcoxon + Cliff's delta 분석 시작...\n")
  
  ## 1) long format + 상대 abundance
  mdf_lvl <- df %>%
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
        TRUE                      ~ NA_character_
      )
    ) %>%
    filter(!is.na(Count), !is.na(Group)) %>%
    group_by(Sample) %>%
    mutate(RelAbund = Count / sum(Count, na.rm = TRUE)) %>%
    ungroup()
  
  cat("  - mdf_lvl rows:", nrow(mdf_lvl), "\n")
  
  ## 2) Wilcoxon + Cliff's delta 계산
  res_delta <- expand.grid(
    Taxon  = unique(mdf_lvl[[lvl]]),
    Group2 = compare_groups,
    stringsAsFactors = FALSE
  ) %>%
    rowwise() %>%
    mutate(
      ref_vals = list(mdf_lvl$RelAbund[mdf_lvl[[lvl]] == Taxon & mdf_lvl$Group == reference_group]),
      cmp_vals = list(mdf_lvl$RelAbund[mdf_lvl[[lvl]] == Taxon & mdf_lvl$Group == Group2]),
      n_ref = length(unlist(ref_vals)),
      n_cmp = length(unlist(cmp_vals)),
      p = if (n_ref >= 2 & n_cmp >= 2) {
        wilcox.test(unlist(ref_vals), unlist(cmp_vals), exact = FALSE)$p.value
      } else NA_real_,
      delta = if (n_ref >= 2 & n_cmp >= 2) {
        effsize::cliff.delta(unlist(cmp_vals), unlist(ref_vals))$estimate
      } else NA_real_,
      mean_abund = mean(c(unlist(ref_vals), unlist(cmp_vals)), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      sig    = !is.na(p) & p < 0.1,
      dir    = ifelse(delta >= 0, "Positive association", "Negative association"),
      size_v = log10(pmax(mean_abund, 1e-6))  # NaN 방지
    )
  
  cat("  - res_delta rows:", nrow(res_delta), "\n")
  cat("  - non-NA p:", sum(!is.na(res_delta$p)), "\n")
  cat("  - raw p < 0.05 개수:", sum(res_delta$p < 0.05, na.rm = TRUE), "\n")
  
  ## 3) 유의한 pairs 필터링
  sig_pairs <- res_delta %>%
    filter(sig, !is.na(delta), !is.na(size_v)) %>%
    select(Taxon, Group2, delta, dir, size_v)
  
  cat("  - sig (Taxon, Group2) pairs:", nrow(sig_pairs), "\n")
  
  if (nrow(sig_pairs) == 0) {
    cat("  >>", lvl, ": 유의한 Taxon-Group2 조합 없음. 스킵.\n")
    next
  }
  
  ## 4) 샘플 단위 데이터 (size_cat 추가)
  plot_df_sample <- sig_pairs %>%
    inner_join(
      mdf_lvl %>%
        filter(Group %in% compare_groups) %>%
        select(Taxon = all_of(lvl), Sample, Group, RelAbund),
      by = c("Taxon", "Group2" = "Group")
    ) %>%
    mutate(
      Sample = factor(Sample, levels = sort(unique(Sample))),
      # size_cat 추가: RelAbund 기준 5단계 분류
      size_cat = cut(
        RelAbund,
        breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 1e-1, Inf),
        labels = c("very low", "low", "mid", "high", "very high"),
        right = FALSE
      )
    )
  
  cat("  - plot_df_sample rows:", nrow(plot_df_sample), "\n")
  
  ## 5) 그룹 평균 데이터
  plot_df_mean <- plot_df_sample %>%
    group_by(Taxon, Group2) %>%
    summarise(
      RelAbund = mean(RelAbund, na.rm = TRUE),
      delta    = first(delta),
      dir      = first(dir),
      .groups  = "drop"
    ) %>%
    mutate(Group2 = factor(Group2, levels = compare_groups))
  
  cat("  - plot_df_mean rows:", nrow(plot_df_mean), "\n")
  
  ## 6) 샘플 피겨
  if (do_sample_plot && nrow(plot_df_sample) > 0) {
    p_tri_sample <- ggplot(plot_df_sample, aes(x = Taxon, y = Sample)) +
      geom_point(
        aes(shape = dir, fill = delta, size = size_cat),
        color = "black", stroke = 0.3, alpha = 0.8
      ) +
      scale_shape_manual(
        values = c("Positive association" = 24, "Negative association" = 25),
        name = "Direction"
      ) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", name = "Effect size") +
      scale_size_manual(
        name = "Relative abundance",
        values = c(6, 10, 14, 18, 22),
        labels = c("very low", "low", "mid", "high", "very high"),
        guide = guide_legend(
          override.aes = list(shape = 24, fill = "black", colour = "black")
        )
      ) +
      theme_bw() +
      theme(
        axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 11),
        axis.title   = element_text(size = 14, face = "bold"),
        plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold")
      ) +
      labs(
        title = paste(lvl, "significant changes vs", reference_group, "(sample-level)"),
        x = lvl, y = "Sample (Direct1~5, Indirect1~5)"
      )
    
    print(p_tri_sample)
    ggsave(file.path(out_dir, paste0("Triangle_SAMPLE_", lvl, "_vs_", reference_group, ".png")),
           p_tri_sample, width = 14, height = 10, dpi = 600)
    write.csv(plot_df_sample, file.path(out_dir, paste0("Triangle_SAMPLE_TABLE_", lvl, "_vs_", reference_group, ".csv")), row.names = FALSE)
    cat("  - 샘플 피겨/CSV 저장 완료\n")
  }
  
  ## 7) 평균 피겨
  if (do_mean_plot && nrow(plot_df_mean) > 0) {
    p_tri_mean <- ggplot(plot_df_mean, aes(x = Taxon, y = Group2)) +
      geom_point(aes(shape = dir, fill = delta, size = RelAbund),
                 color = "black", stroke = 0.7) +
      scale_shape_manual(values = c("Positive association" = 24, "Negative association" = 25),
                         name = "Direction") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                           name = "Effect size (Cliff's delta)") +
      scale_size(range = c(3, 7), name = "Mean relative abundance") +
      theme_bw() +
      theme(
        axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 11),
        axis.title   = element_text(size = 14, face = "bold"),
        plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold")
      ) +
      labs(title = paste(lvl, "significant changes vs", reference_group, "(group mean)"),
           x = lvl, y = "Group")
    
    print(p_tri_mean)
    ggsave(file.path(out_dir, paste0("Triangle_MEAN_", lvl, "_vs_", reference_group, ".png")),
           p_tri_mean, width = 14, height = 10, dpi = 600)
    write.csv(plot_df_mean, file.path(out_dir, paste0("Triangle_MEAN_TABLE_", lvl, "_vs_", reference_group, ".csv")), row.names = FALSE)
    cat("  - 평균 피겨/CSV 저장 완료\n")
  }
}

cat("\n[TRIANGLE] 모든 level 처리 완료.\n")