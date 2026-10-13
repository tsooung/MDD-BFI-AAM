metadata_group <- "metadata1.txt"
level_s <- c("HC", "MDD")
pl_list <- c("HC  MDD")
comparlist <- list(c("HC", "MDD"))

pl2 = c("HC", "MDD")
cl = "species"
pl = "HC  MDD.hebing"
index_list = c("Shannon", "Simpson", "InvSimpson")

library(ggplot2)
library(ggsci)
library(ggsignif)

for (index1 in index_list) {
  metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  metadata00 <- subset(metadata, Group %in% pl2)
  metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
  otutab = read.table(paste0("tem/alpha.", cl, ".", pl, ".csv"), header = T, row.names = 1, sep = ",", quot = "", comment.char = "")

  otutab <- otutab[rownames(metadata00), ]

  df = cbind(otutab[rownames(metadata00), index1], metadata00)
  colnames(df) = c("id", "group", "gender", "age")

  newdir15 = paste0("custom_plots/")
  dir.create(newdir15, recursive = TRUE)

  df$group <- factor(df$group,
                     levels = level_s, ordered = TRUE)

  fig <- ggplot(df, aes(x = group, y = id, color = group)) +
    geom_jitter(aes(color = group), position = position_jitter(0.2), size = 5) +
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black") +
    stat_summary(fun = mean, geom = "point", color = "black", size = 1.5, aes(group = group)) +
    geom_signif(
      comparisons = list(c("HC", "MDD")),
      test = "wilcox.test",
      map_signif_level = TRUE,
      color = "black"
    ) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(colour = NULL) + labs(x = "") +
    ylab(paste0(index1, " index in Species level")) +
    theme(axis.title = element_text(size = 12)) +
    scale_color_aaas(a = 0.55) +
    scale_fill_aaas(a = 0.55) +
    theme_bw() +
    theme(strip.text.y = element_text(face = "italic")) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
    theme(axis.title = element_text(face = "bold", size = 15),
          legend.title = element_text(face = "bold"),
          axis.text.y = element_text(family = "sans", size = 12, colour = 'black'),
          axis.text.x = element_text(size = 12, family = "sans", colour = 'black'),
          legend.text = element_text(family = "sans", size = 12, colour = 'black'),
          legend.text.align = 0)

  fig

  ggsave(path = newdir15, filename = paste0(index1, ".pdf"), fig, width = 120, height = 120, units = "mm")
  ggsave(path = newdir15, filename = paste0(index1, ".png"), fig, width = 120, height = 120, units = "mm")
}

library(vegan)
library(ggsci)
library(ggplot2)
library(ggExtra)

cl_list = c("species")
for (cl in cl_list) {
  print(cl)
  for (pl in pl_list) {
    print(pl)
    pl2 <- as.character(unlist(strsplit(pl, split = "  ")))
    if (length(pl2) != 2) {
      metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
      otutab = read.table(paste0("merged_filtered/out_", cl, "_x100_0.01_20filtered.txt"), header = T, row.names = 1, sep = "\t", quot = "", comment.char = "", check.names = F)

      source("scripts/pcoa_pic.R")

      metadata00 <- subset(metadata, Group %in% pl2)
      metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
      otutab <- otutab[, rownames(metadata00)]

      newdir11 = paste0("custom_plots")
      dir.create(newdir11, recursive = TRUE)

      otutab = otutab
      metadata = metadata00
      cl = cl
      Group = "Group"
      box = "TRUE"
      method = "bray"
      color = "aaas"
      label = "F"
      level_s = pl2
      comparlist = comparlist
      newdir = newdir11
      xline = 1
      yline = 2

      otutab <- otutab[, rownames(metadata)]
      otutab <- otutab[which(rowSums(otutab) > 0), ]
      Transp_otu <- t(otutab)

      bray_dis <- vegdist(Transp_otu, method = method)
      bray_dis
      bray_dis1 <- as.matrix(bray_dis)

      library(devtools)
      library(htmlwidgets)
      library(tidyverse)
      library(plotly)

      p_list = c("ggplot2", "vegan", "ggrepel")
      for (p in p_list) {
        library(p, character.only = T)
      }

      idx = rownames(metadata) %in% rownames(bray_dis1)

      metadata = metadata[idx, , drop = F]
      metadata1 = as.data.frame(lapply(metadata, as.character), stringsAsFactors = F)
      rownames(metadata1) <- rownames(metadata)

      bray_dis1 = bray_dis1[rownames(metadata1), rownames(metadata1)]

      sampFile = as.data.frame(metadata[, Group], row.names = row.names(metadata1))
      sampFile

      pcoa = cmdscale(bray_dis1, k = 3, eig = T)

      points = as.data.frame(pcoa$points)
      points
      eig = pcoa$eig

      points1 = cbind(points, sampFile)
      points1
      colnames(points1) = c("V1", "V2", "V3", Group)
      points1$Group <- factor(points1$Group,
                              levels = level_s, ordered = TRUE)

      otu.adonis <- adonis2(t(otutab) ~ Group, data = metadata, permutations = 999, method = "bray")

      p0 <- ggplot(points1, aes(x = V1, y = V2, color = Group)) +
        geom_point(alpha = 0.68, size = 3, shape = 16) +
        labs(x = paste("PCo ", xline, " (", format(100 * eig[xline] / sum(eig), digits = 4), "%)", sep = ""),
             y = paste("PCo ", yline, " (", format(100 * eig[yline] / sum(eig), digits = 4), "%)", sep = ""), color = Group) +
        theme_bw() +
        theme(
          plot.title = element_text(size = 20, face = "bold", colour = 'black'),
          axis.text.y = element_text(family = "sans", size = 12, colour = 'black'),
          axis.text.x = element_text(size = 12, family = "sans", colour = 'black'),
          axis.title = element_text(size = 12, face = "bold", colour = 'black'),
          legend.text = element_text(family = "sans", size = 9, colour = 'black'),
          legend.position = 'left') +
        stat_ellipse(level = 0.10) +
        scale_color_aaas(a = 0.55) +
        scale_fill_aaas(a = 0.55)
      (p1 <- p0 + ggtext::geom_richtext(
        aes(vjust = 1.1, hjust = 1.05,
            x = Inf, y = Inf,
            label = paste('PERMANOVA: <p>\n  p-value = ',
                          otu.adonis$`Pr(>F)`[1], sep = "")),
        size = 3, family = "sans", colour = 'black'))

      p7 <- ggMarginal(p1, type = "densigram", groupColour = TRUE, groupFill = TRUE, alpha = 0.5)
      ggsave(path = "custom_plots", filename = paste0("PCoA of ", cl, ".7.pdf"), p7, width = 160, height = 120, units = "mm")
      ggsave(path = "custom_plots", filename = paste0("PCoA of ", cl, ".7.png"), p7, width = 160, height = 120, units = "mm")
    }
  }
}

library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer)

cl <- "species"
pl <- "HC  MDD"
pl2 <- as.character(unlist(strsplit(pl, split = "  ")))

path_otu <- paste0("merged_filtered/out_", cl, "_x100_0.01_20filtered.txt")
path_meta <- "metadata1.txt"
path_anno <- "taxnomy/taxnomy_zong.txt"

tax_col <- c("d__Bacteria" = "#F8C4B4",
             "d__Eukaryota" = "#FFEDD5",
             "d__Archaea" = "#CDE990",
             "d__Viruses" = "#99B8B1")

outdata <- read.table(path_otu, header = T, row.names = 1, sep = "\t", check.names = F)
metadata00 <- read.table(path_meta, header = T, row.names = 1, sep = "\t", check.names = F)
anno <- read.table(path_anno, header = F, sep = "|", check.names = F)
names(anno)[1] <- "Kingdom"

anno$Feature <- gsub("^s__", "", anno$V8)
anno1 <- anno[, c("Feature", "Kingdom")]

anno1 <- anno1[anno1$Feature %in% rownames(outdata), ]

common_samples <- intersect(colnames(outdata), rownames(metadata00))
outdata00 <- outdata[, common_samples]
metadata00 <- metadata00[common_samples, ]

available_kingdoms <- intersect(unique(anno1$Kingdom), names(tax_col))

abund_cutoff <- 0.01
prev_cutoff  <- 0.2

lt_list2 <- list()
lt2 <- list()

for (k in available_kingdoms) {
  lt_list2[[k]] <- list()
  k_species <- anno1$Feature[anno1$Kingdom == k]

  for (group_name in pl2) {
    group_samples <- rownames(metadata00[metadata00$Group == group_name, ])

    sub_out <- outdata00[intersect(rownames(outdata00), k_species), group_samples, drop = F]

    if (ncol(sub_out) > 0) {
      mean_abund <- rowMeans(sub_out)
      prevalence <- rowSums(sub_out > 0) / ncol(sub_out)

      keep_idx <- which(mean_abund >= abund_cutoff & prevalence >= prev_cutoff)
      present_species <- names(keep_idx)
    } else {
      present_species <- character(0)
    }

    lt_list2[[k]][[group_name]] <- present_species

    lt2[[group_name]] <- unique(c(lt2[[group_name]], present_species))
  }
}

length(union(lt2$HC, lt2$MDD))

m2 <- make_comb_mat(lt2)

m_list2 <- lapply(lt_list2, make_comb_mat)

rowAnn_list <- list()
for (k in names(tax_col)) {
  if (k %in% names(m_list2)) {
    rowAnn_list[[k]] <- set_size(m_list2[[k]])
  } else {
    rowAnn_list[[k]] <- rep(0, length(pl2))
  }
}
rowAnn_mt2 <- as.matrix(do.call(cbind, rowAnn_list))

colAnn_mt2 <- matrix(0, nrow = length(comb_name(m2)), ncol = length(tax_col),
                     dimnames = list(comb_name(m2), names(tax_col)))

for (k in names(m_list2)) {
  sizes <- comb_size(m_list2[[k]])
  common_combs <- intersect(names(sizes), rownames(colAnn_mt2))
  colAnn_mt2[common_combs, k] <- sizes[common_combs]
}

upset_plot <- UpSet(
  m2,
  set_order = pl2,
  comb_order = order(comb_size(m2), decreasing = TRUE),
  pt_size = unit(5, "mm"),
  lwd = 3,
  bg_col = c("#F0F0F0", "#FFFFFF"),

  top_annotation = HeatmapAnnotation(
    "Intersection size" = anno_barplot(
      colAnn_mt2,
      gp = gpar(fill = tax_col, col = NA),
      height = unit(4, "cm"),
      border = FALSE
    ),
    annotation_name_side = "left",
    annotation_name_rot = 90
  ),

  right_annotation = rowAnnotation(
    "Set size" = anno_barplot(
      rowAnn_mt2,
      gp = gpar(fill = tax_col, col = NA),
      width = unit(3, "cm"),
      border = FALSE
    )
  ),
  row_names_gp = gpar(fontface = "bold", fontsize = 12)
)

draw(upset_plot)

save_dir = paste0("species_analysis/species_composition/upset/")
pdf_path <- paste0(save_dir, "UpSet_Complex_Taxa.pdf")

pdf(file = pdf_path, width = 5, height = 4.5)

draw(upset_plot)

dev.off()

png_path <- paste0(save_dir, "UpSet_Complex_Taxa.png")
png(file = png_path, width = 1200, height = 1000, res = 300)
draw(upset_plot)
dev.off()
