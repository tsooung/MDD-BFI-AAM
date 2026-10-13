print("ReporterScore")

metadata_group <- "metadata1.txt"
pl_list <- c("HC  MDD")

library(ggplot2)
library("ReporterScore")

for (pl in pl_list) {
  pl2 <- as.character(unlist(strsplit(pl, split = "  ")))

  metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  metadata00 <- subset(metadata, Group %in% pl2)

  KO_abun <- read.table("functional_maaslin/out_KO.txt", header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  KO_abun <- KO_abun[, colnames(KO_abun) %in% row.names(metadata00)]
  KO_abun <- KO_abun[which(rowSums(KO_abun) > 0), ]
  metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]

  newdir27 = paste0("functional_maaslin/reporter_score/", pl)
  dir.create(newdir27, recursive = TRUE)

  if (length(pl2) == 2) {
    reporter_res <- GRSA(KO_abun, "Group", metadata00,
                         mode = "directed",
                         method = "wilcox.test", perm = 999,
                         type = "pathway", feature = "ko"
    )
  }
  ko_stat <- reporter_res$ko_stat
  write.table(ko_stat, file = paste0(newdir27, "/ko_stat.txt"), sep = "\t", quote = F, row.names = F)

  reporter_s <- reporter_res$reporter_s
  write.table(reporter_s, file = paste0(newdir27, "/reporter_s.txt"), sep = "\t", quote = F, row.names = F)

  modulelist <- reporter_res$modulelist
  write.table(modulelist, file = paste0(newdir27, "/pathwaylist.txt"), sep = "\t", quote = F, row.names = F)

  (p1 <- plot_report_bar(reporter_res, rs_threshold = c(-2.5, 2.5), facet_level = TRUE))
  ggsave(path = newdir27, filename = paste0("report_bar.pdf"), p1, width = 270, height = 170, units = "mm")
  ggsave(path = newdir27, filename = paste0("report_bar.png"), p1, width = 270, height = 170, units = "mm")

}

library(tidyverse)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggsci)

df <- read_tsv("functional_maaslin/maaslin2_ec/ec.txt", col_types = cols(
  id     = col_character(),
  module = col_character(),
  ecid   = col_character(),
  anno   = col_character(),
  pval   = col_double(),
  coef   = col_double()
))

df <- df %>%
  mutate(
    log10p = -log10(pval),
    anno_short = str_replace(anno, "\\s*\\([^)]*\\)", ""),
    sig = case_when(
      pval < 0.001 ~ "***",
      pval < 0.01  ~ "**",
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    )
  ) %>%
  arrange(desc(log10p))

df_plot <- df %>%
  mutate(
    neg_log_pval = -log10(pval),
    y_label = paste0(ecid, "  |  ", str_trunc(anno, 60, "right", "..."))
  ) %>%
  group_by(module) %>%
  arrange(desc(neg_log_pval), .by_group = TRUE) %>%
  mutate(y_label = factor(y_label, levels = rev(y_label))) %>%
  ungroup()

p1 <- ggplot(data = df_plot,
             aes(x = neg_log_pval, y = y_label, fill = module)) +

  geom_bar(stat = "identity", width = 0.7) +

  scale_fill_aaas(alpha = 0.6) +

  labs(
    x = "-log10(p-value)",
    y = "",
    title = "EC Enrichment"
  ) +

  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 9, family = "sans"),
    axis.text.y = element_text(size = 10, family = "sans", color = "black"),

    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),

    text = element_text(family = "sans", size = 10),

    panel.background = element_rect(fill = "grey92"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(
      colour = "grey50",
      size = 0.2,
      linetype = 2,
      lineend = 1
    ),

    legend.key = element_rect(fill = "white"),

    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(family = "sans", size = 11, face = "bold")
  )

p1

ggsave(paste0("EC.diff.1.pdf"), p1, width = 200, height = 150, units = "mm")
ggsave(paste0("EC.diff.1.png"), p1, width = 200, height = 150, units = "mm")

library(readxl)
library(openxlsx)
library(tidyverse)
library(ppcor)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

pheno <- read_xlsx("phenotype_mantel/phenotype_data.xlsx", sheet = "Sheet3") %>%
  column_to_rownames("ID")

module_diff <- read_xlsx(
  "functional_maaslin/grsa_enrichment/grsa_module_pathway.xlsx",
  sheet = "raw_results"
) %>% column_to_rownames("ID")

pathway_filter <- read.table(
  "functional_maaslin/out_Module.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  .[rownames(module_diff), ]

metadata <- read.table(
  "metadata1.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  stringsAsFactors = FALSE
)

pheno <- pheno[rownames(metadata), ]
pathway_filter <- t(pathway_filter[rowSums(pathway_filter) > 0, ])
pathway_filter <- pathway_filter[rownames(pheno), ]

covariates <- metadata[, c("gender", "ages", "BMI")]
covariates$gender <- as.numeric(factor(covariates$gender))

cor_res <- data.frame()

for (mic in colnames(pheno)) {
  for (met in colnames(pathway_filter)) {

    x <- as.numeric(pheno[, mic])
    y <- as.numeric(pathway_filter[, met])

    idx <- complete.cases(x, y, covariates)

    if (sum(idx) <= ncol(covariates) + 2) next

    tmp <- pcor.test(
      x[idx],
      y[idx],
      covariates[idx, ],
      method = "spearman"
    )

    cor_res <- rbind(
      cor_res,
      data.frame(
        mic = mic,
        met = met,
        estimate = tmp$estimate,
        p.value = tmp$p.value
      )
    )
  }
}

dir.create("module_phenotype_correlation", showWarnings = FALSE)

write.xlsx(cor_res, paste0("module_phenotype_correlation/correlation_unfiltered.all.xlsx"),
           append = FALSE,
           col.names = TRUE,
           row.names = FALSE)

cor_res <- read.xlsx("module_phenotype_correlation/correlation_unfiltered.all.xlsx", sheetIndex = 1)
cor_res$qval <- p.adjust(cor_res$p.value, method = "BH")
write.xlsx(cor_res, paste0("module_phenotype_correlation/correlation_unfiltered.all_addqvalue.xlsx"),
           append = FALSE,
           col.names = TRUE,
           row.names = FALSE)

plotmat <- cor_res %>%
  dplyr::select(mic, met, estimate) %>%
  tidyr::pivot_wider(names_from = met, values_from = estimate) %>%
  column_to_rownames("mic") %>%
  as.matrix()

stars <- cor_res %>%
  dplyr::select(mic, met, qval) %>%
  tidyr::pivot_wider(names_from = met, values_from = qval) %>%
  column_to_rownames("mic") %>%
  as.matrix()

stars_char <- matrix("", nrow = nrow(plotmat), ncol = ncol(plotmat))
rownames(stars_char) <- rownames(plotmat)
colnames(stars_char) <- colnames(plotmat)

for (i in 1:nrow(plotmat)) {
  for (j in 1:ncol(plotmat)) {
    if (abs(plotmat[i, j]) >= 0) {
      if (stars[i, j] < 0.01) {
        stars_char[i, j] <- "**"
      } else if (stars[i, j] < 0.05) {
        stars_char[i, j] <- "*"
      }
    }
  }
}

anno_df <- module_diff[
  colnames(plotmat),
  c("ReporterScore", "Significant_up_num")
]

anno_df$Enrichment <- ifelse(
  anno_df$ReporterScore > 0,
  "MDD",
  "HC"
)

writeLines(colnames(plotmat)[column_order(p)],
           "module_phenotype_correlation/module_column_order_padj.txt")
writeLines(rownames(plotmat)[row_order(p)],
           "module_phenotype_correlation/module_row_order_padj.txt")

plotmat_t <- t(plotmat)

stars_char_t <- t(stars_char)

anno_df <- anno_df[order(rownames(anno_df)), ]
plotmat_t <- plotmat_t[rownames(plotmat_t) %in% rownames(anno_df), ]
rownames(plotmat_t) <- rownames(anno_df)

max_abs <- max(abs(anno_df$ReporterScore), na.rm = TRUE)

reporter_col_fun <- colorRamp2(
  c(-max_abs, 0, max_abs),
  c("#31688E", "white", "#FB9A99")
)

cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(
    stars_char_t[i, j],
    x = x,
    y = y - height * 0.3,
    gp = gpar(fontsize = 12, col = "black")
  )
}

library(grid)

anno_df$ReporterScore_abs <- abs(anno_df$ReporterScore)

max_abs <- max(abs(anno_df$ReporterScore), na.rm = TRUE)
reporter_col_fun <- colorRamp2(
  c(-max_abs, 0, max_abs),
  c("#1B9E77", "white", "#D95F02")
)

right_bar <- rowAnnotation(
  ReporterScore = anno_barplot(
    anno_df$ReporterScore_abs,
    baseline = 0,
    bar_width = 0.8,
    gp = gpar(
      fill = reporter_col_fun(anno_df$ReporterScore),
      col = NA
    )
  ),
  width = unit(1.5, "cm"),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 10),
  border = FALSE
)

reporter_legend <- Legend(
  title = "ReporterScore",
  col_fun = reporter_col_fun,
  at = c(-max_abs, -max_abs / 2, 0, max_abs / 2, max_abs),
  labels = c(
    sprintf("-%.1f", max_abs),
    sprintf("-%.1f", max_abs / 2),
    "0",
    sprintf("%.1f", max_abs / 2),
    sprintf("%.1f", max_abs)
  ),
  direction = "vertical",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9),
  grid_height = unit(0.8, "cm"),
  grid_width = unit(0.6, "cm"),
  border = "transparent"
)

p_transposed <- Heatmap(
  plotmat_t,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  right_annotation = right_bar,
  cell_fun = cell_fun,
  rect_gp = gpar(col = "grey", lwd = 0.5),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)

pdf("module_phenotype_correlation/pathway_pheno_final_padj.pdf",
    width = 8,
    height = 6)

draw(
  p_transposed,
  annotation_legend_list = list(reporter_legend),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

png("module_phenotype_correlation/pathway_pheno_final_padj.png",
    width = 2400,
    height = 1800,
    res = 300,
    bg = "white")

draw(
  p_transposed,
  annotation_legend_list = list(reporter_legend),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()
