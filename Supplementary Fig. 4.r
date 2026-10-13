metadata_group <- "metadata_paired.txt"
level_s <- c("MDD-Pre", "MDD-Post")
pl_list <- c("MDD-Pre  MDD-Post")
comparlist <- list(c("MDD-Pre", "MDD-Post"))

cl <- "species"
pl <- "MDD-Pre  MDD-Post"
pl2 <- c("MDD-Pre", "MDD-Post")
newdir15 <- "custom_plots/"

library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(here)
library(PupillometryR)

dir.create(newdir15, recursive = TRUE)

index_list = c("Shannon", "Simpson", "InvSimpson")
for (index1 in index_list) {
  metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  metadata00 <- subset(metadata, Group %in% pl2)
  metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
  otutab = read.table(paste0("tem/alpha.", cl, ".", pl, ".csv"), header = T, row.names = 1, sep = ",", quot = "", comment.char = "")

  otutab <- otutab[rownames(metadata00), ]

  df = cbind(otutab[rownames(metadata00), index1], metadata00)
  colnames(df) = c("id", "group", "gender", "age")

  df$group <- factor(df$group,
                     levels = level_s, ordered = TRUE)

  print(index1)

  p <- ggplot(df, aes(group, id, fill = group)) +
    geom_flat_violin(aes(fill = group), position = position_nudge(x = .25), color = "black", alpha = 0.4) +
    geom_jitter(position = position_jitter(0.15), size = 4, alpha = 0.4, aes(color = group, shape = group), na.rm = TRUE, inherit.aes = TRUE) +
    geom_boxplot(width = .1, position = position_nudge(x = 0.25), fill = "white", size = 0.5) +
    geom_signif(
      comparisons = list(c("MDD-Pre", "MDD-Post")),
      test = "wilcox.test",
      map_signif_level = TRUE,
      color = "black"
    ) +
    theme(axis.text = element_text(size = 12),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(family = "sans", colour = 'black', size = 15),
          axis.line = element_line(colour = 'black', size = 0.75),
          legend.position = "top",
          text = element_text(family = "sans", size = 10),
          panel.background = element_rect(fill = 'white'),
          legend.key = element_rect(fill = 'white')) +

    labs(x = "", y = paste0(index1, " index in Species level")) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 20, face = "bold", colour = 'black'),
      legend.position = "bo",
      panel.grid = element_blank(),
      text = element_text(family = "sans", size = 10, colour = 'black'),
      axis.text.y = element_text(family = "sans", size = 20, colour = 'black'),
      axis.text.x = element_text(size = 20, hjust = 0.5, colour = "black",
                                 vjust = 0.5, family = "sans"))

  p

  p = p +
    scale_color_aaas(a = 0.55) +
    scale_fill_aaas(a = 0.55)

  ggsave(path = newdir15, filename = paste0(index1, ".1.pdf"), p, width = 130, height = 140, units = "mm")
  ggsave(path = newdir15, filename = paste0(index1, ".1.png"), p, width = 130, height = 140, units = "mm")

  p

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

      p5 <- ggMarginal(p1, type = "boxplot", groupColour = TRUE, groupFill = TRUE, alpha = 0.5)
      ggsave(path = "custom_plots", filename = paste0("PCoA of ", cl, ".5.pdf"), p5, width = 160, height = 120, units = "mm")
      ggsave(path = "custom_plots", filename = paste0("PCoA of ", cl, ".5.png"), p5, width = 160, height = 120, units = "mm")
    }
  }
}
