library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(stringr)

micro_file <- "User microbiome data.csv"
fa_file <- "User fatty acid data.csv"
out_dir <- "User workstation pathway"

micro <- read.csv(micro_file, check.names = FALSE)
fa <- read.csv(fa_file, row.names = 1, check.names = FALSE)

clean_ids <- function(x) {
  x %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("[\\._\\-]", "") %>%
    str_replace_all("Sample", "Sample") %>%  #Target name (direct or indirect any other)
    str_replace_all("ctrl", "Sample")
}

rownames(fa) <- clean_ids(rownames(fa))

target_samples <- c("Sample 1", " Sample 1", " Sample 1")   #Target name (direct or indirect any other)

stopifnot(all(target_samples %in% rownames(fa)))

fa_sel <- fa[target_samples, , drop = FALSE]
fa_mat <- t(fa_sel)

stopifnot(all(c("FullTaxonomy", "Name") %in% colnames(micro)))

rownames(micro) <- micro$Name

micro <- micro %>%
  dplyr::select(-FullTaxonomy, -Name)

stopifnot(all(target_samples %in% colnames(micro)))

micro_mat <- micro[, target_samples, drop = FALSE]

common_samples <- intersect(
  colnames(micro_mat),
  colnames(fa_mat)
)

stopifnot(length(common_samples) >= 2)

micro_mat <- micro_mat[, common_samples, drop = FALSE]
fa_mat <- fa_mat[, common_samples, drop = FALSE]

micro_mat_before <- micro_mat
micro_mat[is.na(micro_mat)] <- 0

nonzero_features <- rowSums(micro_mat) > 0
variable_features <- apply(micro_mat, 1, sd) > 0

micro_mat <- micro_mat[
  nonzero_features & variable_features,
  ,
  drop = FALSE
]

removed_features <- setdiff(
  rownames(micro_mat_before),
  rownames(micro_mat)
)

write.csv(
  removed_features,
  file.path(
    out_dir,
    "Removed_Microbiome_Features_Direct.csv"
  ),
  row.names = FALSE
)

cor_res <- expand.grid(
  micro = rownames(micro_mat),
  fa = rownames(fa_mat),
  stringsAsFactors = FALSE
) %>%
  mutate(
    r = NA_real_,
    p = NA_real_
  )

for (i in seq_len(nrow(cor_res))) {
  g <- cor_res$micro[i]
  f <- cor_res$fa[i]

  x <- as.numeric(micro_mat[g, ])
  y <- as.numeric(fa_mat[f, ])

  if (sd(x) > 0 && sd(y) > 0) {
    test <- suppressWarnings(
      cor.test(
        x,
        y,
        method = "pearson"
      )
    )

    cor_res$r[i] <- unname(test$estimate)
    cor_res$p[i] <- test$p.value
  }
}

cor_mat <- cor_res %>%
  dplyr::select(micro, fa, r) %>%
  pivot_wider(
    names_from = fa,
    values_from = r
  ) %>%
  column_to_rownames("micro") %>%
  as.matrix()

p_mat <- cor_res %>%
  dplyr::select(micro, fa, p) %>%
  pivot_wider(
    names_from = fa,
    values_from = p
  ) %>%
  column_to_rownames("micro") %>%
  as.matrix()

cor_mat[!is.finite(cor_mat)] <- 0
p_mat[!is.finite(p_mat)] <- 1

cut_p <- 0.05
sig_matrix <- p_mat < cut_p

row_keep <- rowSums(
  sig_matrix,
  na.rm = TRUE
) > 0

col_keep <- colSums(
  sig_matrix,
  na.rm = TRUE
) > 0

cor_mat_sig <- cor_mat[
  row_keep,
  col_keep,
  drop = FALSE
]

p_mat_sig <- p_mat[
  row_keep,
  col_keep,
  drop = FALSE
]

fatty_acid_order <- c(
  "C12:0",
  "C14:0",
  "C15:0",
  "C16:0",
  "C17:0",
  "C18:0",
  "C20:0",
  "C22:0",
  "C18:2n6t",
  "C18:3n3"
)

keep_cols <- fatty_acid_order[
  fatty_acid_order %in% colnames(cor_mat_sig)
]

cor_mat_sig <- cor_mat_sig[
  ,
  keep_cols,
  drop = FALSE
]

p_mat_sig <- p_mat_sig[
  ,
  keep_cols,
  drop = FALSE
]

rows_no_na <- apply(
  cor_mat_sig,
  1,
  function(x) all(!is.na(x))
)

row_dend <- if (sum(rows_no_na) > 1) {
  as.dendrogram(
    hclust(
      dist(
        cor_mat_sig[
          rows_no_na,
          ,
          drop = FALSE
        ]
      )
    )
  )
} else {
  FALSE
}

p_star_fun <- function(pval) {
  if (is.na(pval)) {
    return("")
  }

  if (pval < 0.001) {
    return("***")
  }

  if (pval < 0.01) {
    return("**")
  }

  if (pval < 0.05) {
    return("*")
  }

  ""
}

col_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#084eff", "white", "#eb1717")
)

ht <- Heatmap(
  cor_mat_sig,
  name = "Pearson (User sample)", #Target name (direct or indirect any other)
  col = col_fun,
  cluster_rows = row_dend,
  cluster_columns = FALSE,
  column_order = keep_cols,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  border = TRUE,
  rect_gp = gpar(
    col = "white",
    lwd = 0.5
  ),
  row_names_gp = gpar(
    fontsize = 7
  ),
  column_names_gp = gpar(
    fontsize = 10
  ),
  heatmap_legend_param = list(
    at = c(-1, -0.5, 0, 0.5, 1)
  ),
  layer_fun = function(j, i, x, y, w, h, fill) {
    p_values <- p_mat_sig[cbind(i, j)]
    
    stars <- vapply(
      p_values,
      p_star_fun,
      character(1)
    )
    
    keep <- stars != ""
    
    if (any(keep)) {
      grid.text(
        stars[keep],
        x[keep],
        y[keep] - unit(1, "mm"),
        gp = gpar(
          fontsize = 10,
          fontface = "bold",
          col = "white"
        )
      )
    }
  }
)

draw(
  ht,
  heatmap_legend_side = "right"
)

pdf(
  file.path(
    out_dir,
    "Microbiome_FA_Correlation_Direct1to3_pval_Final.pdf"
  ),
  width = 8,
  height = 6
)

draw(
  ht,
  heatmap_legend_side = "right"
)

dev.off()

png(
  file.path(
    out_dir,
    "Microbiome_FA_Correlation_Direct1to3_pval_Final.png"
  ),
  width = 2000,
  height = 2700,
  res = 300
)

draw(
  ht,
  heatmap_legend_side = "right"
)

dev.off()

cor_table_sig <- cor_res %>%
  filter(
    micro %in% rownames(cor_mat_sig),
    fa %in% colnames(cor_mat_sig)
  ) %>%
  mutate(
    Significant = !is.na(p) & p < cut_p
  ) %>%
  arrange(p) %>%
  dplyr::select(
    Microbe = micro,
    FattyAcid = fa,
    Pearson_r = r,
    P_value = p,
    Significant
  )

write.csv(
  cor_table_sig,
  file.path(
    out_dir,
    "Microbiome_FA_Correlation_select direct or indirect_Table_Final.csv"
  ),
  row.names = FALSE
)