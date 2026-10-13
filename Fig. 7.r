rm(list = ls())

metadata_group <- "metadata_paired.txt"
level_s <- c("MDD-Pre", "MDD-Post")
pl_list <- c("MDD-Pre  MDD-Post")
comparlist <- list(c("MDD-Pre", "MDD-Post"))
continuous <- c("ages", "BMI")
discrete <- "gender"

dir.create("tem", recursive = TRUE)
dir.create("tex", recursive = TRUE)

print("alpha diversity index calculation")
library(vegan)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(xlsx)

newdir3 = "species_analysis/alpha"
dir.create(newdir3, recursive = TRUE)

for (pl in pl_list) {
  print(pl)

  newdir3 = paste0("species_analysis/alpha/", pl)
  dir.create(newdir3, recursive = TRUE)

  fit <- try(file.remove(paste0(newdir3, "/alpha.p.xlsx")))
  if ("try-error" %in% class(fit)) {
    next
  } else {
    file.remove(paste0(newdir3, "/alpha.p.xlsx"))
    file.remove(paste0(newdir3, "/alpha.xlsx"))
  }

  cl_list = c("species", "genus", "class", "family", "order", "phylum")
  for (cl in cl_list) {
    print(cl)
    pl2 <- as.character(unlist(strsplit(pl, split = "  ")))
    metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
    metadata00 <- subset(metadata, Group %in% pl2)
    metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
    otu_table <- read.delim(paste0("kraken2/out_kraken2_gtdb_bacteria_", cl, ".txt"), header = T, sep = "\t", row.names = 1, quot = "", stringsAsFactors = FALSE, check.names = F)
    otu_table <- otu_table[, rownames(metadata00)]
    if (length(pl2) == 2) {
      newdir3 = paste0("species_analysis/alpha/", pl)
      dir.create(newdir3, recursive = TRUE)
      otu = t(otu_table)
      Shannon = vegan::diversity(otu, index = "shannon")
      Simpson = vegan::diversity(otu, index = "simpson")
      InvSimpson = vegan::diversity(otu, index = "invsimpson")
      alpha_div = cbind(Shannon, Simpson, InvSimpson)
      write.csv("Sample\t", file = paste0("tem/alpha.", cl, ".", pl, ".csv"), append = F, quote = F, eol = "", row.names = F, col.names = F)
      suppressWarnings(write.csv(alpha_div, file = paste0("tem/alpha.", cl, ".", pl, ".csv"), append = T, quote = F, sep = ",", eol = "\n", na = "NA", dec = "."))
      library(xlsx)
      write.xlsx(alpha_div, paste0(newdir3, "/alpha.xlsx"),
                 sheetName = cl, append = TRUE,
                 col.names = TRUE,
                 row.names = TRUE)

      alpha_div00 <- merge(alpha_div, metadata00, by = "row.names")

      Shannon_p <- wilcox.test(Shannon ~ Group, data = alpha_div00)
      Simpson_p <- wilcox.test(Simpson ~ Group, data = alpha_div00)
      InvSimpson_p <- wilcox.test(InvSimpson ~ Group, data = alpha_div00)

      diversity_p = cbind(Shannon_p, Simpson_p, InvSimpson_p)
      diversity_p_t <- as.data.frame(t(diversity_p))[, c("statistic", "p.value", "method")]

      write.xlsx(diversity_p_t, paste0(newdir3, "/alpha.p.xlsx"),
                 sheetName = cl, append = TRUE,
                 col.names = TRUE,
                 row.names = TRUE)
    }

    if (length(pl2) != 2) {
      otu = t(otu_table)
      Shannon = vegan::diversity(otu, index = "shannon")
      Simpson = vegan::diversity(otu, index = "simpson")
      InvSimpson = vegan::diversity(otu, index = "invsimpson")
      alpha_div = cbind(Shannon, Simpson, InvSimpson)
      write.csv("Sample\t", file = paste0("tem/alpha.", cl, ".", pl, ".csv"), append = F, quote = F, eol = "", row.names = F, col.names = F)
      suppressWarnings(write.csv(alpha_div, file = paste0("tem/alpha.", cl, ".", pl, ".csv"), append = T, quote = F, sep = ",", eol = "\n", na = "NA", dec = "."))

      write.xlsx(alpha_div, paste0(newdir3, "/alpha.xlsx"),
                 sheetName = cl, append = TRUE,
                 col.names = TRUE,
                 row.names = TRUE)

      alpha_div00 <- merge(alpha_div, metadata, by = "row.names")

      Shannon_p <- kruskal.test(Shannon ~ Group, data = alpha_div00)
      Simpson_p <- kruskal.test(Simpson ~ Group, data = alpha_div00)
      InvSimpson_p <- kruskal.test(InvSimpson ~ Group, data = alpha_div00)

      diversity_p = cbind(Shannon_p, Simpson_p, InvSimpson_p)
      diversity_p_t <- as.data.frame(t(diversity_p))[, c("statistic", "p.value", "method")]

      write.xlsx(diversity_p_t, paste0(newdir3, "/alpha.p.xlsx"),
                 sheetName = cl, append = TRUE,
                 col.names = TRUE,
                 row.names = TRUE)
    }
  }
}

print("alpha diversity plotting")

for (pl in pl_list) {
  cl_list = c("species", "genus", "class", "family", "order", "phylum")
  for (cl in cl_list) {
    print(cl)
    p_list = c("Shannon", "Simpson", "InvSimpson")
    for (i in p_list) {
      print(i)

      source("scripts/boxplot_continuous.R")
      pl2 <- as.character(unlist(strsplit(pl, split = "  ")))
      metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
      metadata00 <- subset(metadata, Group %in% pl2)
      metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
      otutab = read.table(paste0("tem/alpha.", cl, ".", pl, ".csv"), header = T, row.names = 1, sep = ",", quot = "", comment.char = "")

      otutab <- otutab[rownames(metadata00), ]

      if (length(pl2) == 2) {
        p = boxplot_continuous(otutab, metadata00, index = i, Group = "Group", color = "NA", level_s = pl2, signif = "T", label = "F", comparlist = list(c(pl2)))
        p1 = boxplot_continuous(otutab, metadata00, index = i, Group = "Group", color = "NA", level_s = pl2, signif = "T", label = "T", comparlist = list(c(pl2)))

        (p1.1 <- p + ylab(paste0(i, " index in ", cl, " level")) + scale_color_aaas(a = 0.55))
        (p2.1 <- p1 + ylab(paste0(i, " index in ", cl, " level")) + scale_color_aaas(a = 0.55))
        newdir4 = paste0("species_analysis/alpha/", pl, "/alpha picture1/", cl)
        dir.create(newdir4, recursive = TRUE)
        ggsave(path = newdir4, filename = paste0(i, " index in ", cl, " level.1.pdf"), p1.1, width = 100, height = 130, units = "mm")
        ggsave(path = newdir4, filename = paste0(i, " index in ", cl, " level.1.png"), p1.1, width = 100, height = 130, units = "mm")
        ggsave(path = newdir4, filename = paste0(i, " index in ", cl, " level.1.point.pdf"), p2.1, width = 100, height = 130, units = "mm")
        ggsave(path = newdir4, filename = paste0(i, " index in ", cl, " level.1.point.png"), p2.1, width = 100, height = 130, units = "mm")

        (p1.2 <- p + ylab(paste0(i, " index in ", cl, " level")) + scale_color_d3(a = 0.55))
        (p2.2 <- p1 + ylab(paste0(i, " index in ", cl, " level")) + scale_color_d3(a = 0.55))
        newdir5 = paste0("species_analysis/alpha/", pl, "/alpha picture2/", cl)
        dir.create(newdir5, recursive = TRUE)
        ggsave(path = newdir5, filename = paste0(i, " index in ", cl, " level.2.pdf"), p1.2, width = 100, height = 130, units = "mm")
        ggsave(path = newdir5, filename = paste0(i, " index in ", cl, " level.2.png"), p1.2, width = 100, height = 130, units = "mm")
        ggsave(path = newdir5, filename = paste0(i, " index in ", cl, " level.2.point.pdf"), p2.2, width = 100, height = 130, units = "mm")
        ggsave(path = newdir5, filename = paste0(i, " index in ", cl, " level.2.point.png"), p2.2, width = 100, height = 130, units = "mm")
      }

      if (length(pl2) != 2) {
        comparlist2 <- list()
        for (o in 1:(length(pl2) - 1)) {
          for (j in (o + 1):length(pl2)) {
            combination <- c(pl2[o], pl2[j])
            comparlist2 <- c(comparlist2, list(combination))
          }
        }
        comparlist2 <- unique(comparlist2)

        p = boxplot_continuous(otutab, metadata00, index = i, Group = "Group", color = "NA", level_s = pl2, signif = "T", label = "F", comparlist = comparlist2)
        p1 = boxplot_continuous(otutab, metadata00, index = i, Group = "Group", color = "NA", level_s = pl2, signif = "T", label = "T", comparlist = comparlist2)

        (p1.1 <- p + ylab(paste0(i, " index in ", cl, " level")) + scale_color_aaas(a = 0.55))
        (p2.1 <- p1 + ylab(paste0(i, " index in ", cl, " level")) + scale_color_aaas(a = 0.55))
        newdir4 = paste0("species_analysis/alpha/", pl, "/alpha picture1/", cl)
        dir.create(newdir4, recursive = TRUE)
        ggsave(path = newdir4, filename = paste0(i, " index in ", cl, " level.1.pdf"), p1.1, width = 100, height = 130, units = "mm")
        ggsave(path = newdir4, filename = paste0(i, " index in ", cl, " level.1.png"), p1.1, width = 100, height = 130, units = "mm")
        ggsave(path = newdir4, filename = paste0(i, " index in ", cl, " level.1.point.pdf"), p2.1, width = 100, height = 130, units = "mm")
        ggsave(path = newdir4, filename = paste0(i, " index in ", cl, " level.1.point.png"), p2.1, width = 100, height = 130, units = "mm")

        (p1.2 <- p + ylab(paste0(i, " index in ", cl, " level")) + scale_color_d3(a = 0.55))
        (p2.2 <- p1 + ylab(paste0(i, " index in ", cl, " level")) + scale_color_d3(a = 0.55))
        newdir5 = paste0("species_analysis/alpha/", pl, "/alpha picture2/", cl)
        dir.create(newdir5, recursive = TRUE)
        ggsave(path = newdir5, filename = paste0(i, " index in ", cl, " level.2.pdf"), p1.2, width = 100, height = 130, units = "mm")
        ggsave(path = newdir5, filename = paste0(i, " index in ", cl, " level.2.png"), p1.2, width = 100, height = 130, units = "mm")
        ggsave(path = newdir5, filename = paste0(i, " index in ", cl, " level.2.point.pdf"), p2.2, width = 100, height = 130, units = "mm")
        ggsave(path = newdir5, filename = paste0(i, " index in ", cl, " level.2.point.png"), p2.2, width = 100, height = 130, units = "mm")
      }
    }
  }
}

cl_list = c("species", "genus", "family", "class", "order", "phylum")
for (cl in cl_list) {
  print(cl)
  for (pl in pl_list) {
    print(pl)
    pl2 <- as.character(unlist(strsplit(pl, split = "  ")))
    if (length(pl2) != 2) {
      metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
      otutab = read.table(paste0("kraken2/out_kraken2_gtdb_bacteria_", cl, ".txt"), header = T, row.names = 1, sep = "\t", quot = "", comment.char = "", check.names = F)
      source("scripts/pcoa_pic.R")

      metadata00 <- subset(metadata, Group %in% pl2)
      metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
      otutab <- otutab[, rownames(metadata00)]

      newdir11 = paste0("species_analysis/beta/pcoa/", pl, "/pcoa picture1")
      dir.create(newdir11, recursive = TRUE)
      p1 = beta_pcoa3(otutab = otutab, metadata00, cl = cl, Group = "Group", box = "FALSE", method = "bray", color = "aaas", label = "F", level_s = pl2, comparlist = comparlist, newdir = newdir11)
      ggsave(path = newdir11, filename = paste0("PCoA of ", cl, ".1.pdf"), p1, width = 150, height = 120, units = "mm")
      ggsave(path = newdir11, filename = paste0("PCoA of ", cl, ".1.png"), p1, width = 150, height = 120, units = "mm")
      p1.2 = beta_pcoa3(otutab = otutab, metadata00, cl = cl, Group = "Group", box = "FALSE", method = "bray", color = "aaas", label = "T", level_s = pl2, comparlist = comparlist, newdir = newdir11)
      ggsave(path = newdir11, filename = paste0("PCoA of ", cl, ".1.point.pdf"), p1.2, width = 150, height = 120, units = "mm")
      ggsave(path = newdir11, filename = paste0("PCoA of ", cl, ".1.point.png"), p1.2, width = 150, height = 120, units = "mm")

      newdir12 = paste0("species_analysis/beta/pcoa/", pl, "/pcoa picture2")
      dir.create(newdir12, recursive = TRUE)
      p2 = beta_pcoa3(otutab = otutab, metadata00, cl = cl, Group = "Group", box = "FALSE", method = "bray", color = "d3", label = "F", level_s = pl2, comparlist = comparlist, newdir = newdir12)
      ggsave(path = newdir12, filename = paste0("PCoA of ", cl, ".2.pdf"), p2, width = 150, height = 120, units = "mm")
      ggsave(path = newdir12, filename = paste0("PCoA of ", cl, ".2.png"), p2, width = 150, height = 120, units = "mm")
      p2.2 = beta_pcoa3(otutab = otutab, metadata00, cl = cl, Group = "Group", box = "FALSE", method = "bray", color = "d3", label = "T", level_s = pl2, comparlist = comparlist, newdir = newdir12)
      ggsave(path = newdir12, filename = paste0("PCoA of ", cl, ".2.point.pdf"), p2.2, width = 150, height = 120, units = "mm")
      ggsave(path = newdir12, filename = paste0("PCoA of ", cl, ".2.point.png"), p2.2, width = 150, height = 120, units = "mm")
    }

    if (length(pl2) == 2) {
      otutab = read.table(paste0("kraken2/out_kraken2_gtdb_bacteria_", cl, ".txt"), header = T, row.names = 1, sep = "\t", comment.char = "", quot = "", check.names = F)
      metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)

      metadata00 <- subset(metadata, Group %in% pl2)
      metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
      otutab <- otutab[, rownames(metadata00)]

      newdir12.1 = paste0("species_analysis/beta/pcoa/", pl, "/pcoa picture1/")
      dir.create(newdir12.1, recursive = TRUE)
      source("scripts/pcoa_pic.R")
      p3 = beta_pcoa3(otutab = otutab, metadata00, cl = cl, Group = "Group", box = "FALSE", method = "bray", color = "aaas", label = "F", level_s = pl2, comparlist = list(c(pl2)), newdir = newdir12.1)
      ggsave(path = newdir12.1, filename = paste0("PCoA of ", cl, ".1.pdf"), p3, width = 150, height = 120, units = "mm")
      ggsave(path = newdir12.1, filename = paste0("PCoA of ", cl, ".1.png"), p3, width = 150, height = 120, units = "mm")
      p3.2 = beta_pcoa3(otutab = otutab, metadata00, cl = cl, Group = "Group", box = "FALSE", method = "bray", color = "aaas", label = "T", level_s = pl2, comparlist = list(c(pl2)), newdir = newdir12.1)
      ggsave(path = newdir12.1, filename = paste0("PCoA of ", cl, ".1.point.pdf"), p3.2, width = 150, height = 120, units = "mm")
      ggsave(path = newdir12.1, filename = paste0("PCoA of ", cl, ".1.point.png"), p3.2, width = 150, height = 120, units = "mm")

      newdir12.2 = paste0("species_analysis/beta/pcoa/", pl, "/pcoa picture2/")
      dir.create(newdir12.2, recursive = TRUE)
      p4 = beta_pcoa3(otutab = otutab, metadata00, cl = cl, Group = "Group", box = "FALSE", method = "bray", color = "d3", label = "F", level_s = pl2, comparlist = list(c(pl2)), newdir = newdir12.2)
      ggsave(path = newdir12.2, filename = paste0("PCoA of ", cl, ".2.pdf"), p4, width = 150, height = 120, units = "mm")
      ggsave(path = newdir12.2, filename = paste0("PCoA of ", cl, ".2.png"), p4, width = 150, height = 120, units = "mm")
      p4.2 = beta_pcoa3(otutab = otutab, metadata00, cl = cl, Group = "Group", box = "FALSE", method = "bray", color = "d3", label = "T", level_s = pl2, comparlist = list(c(pl2)), newdir = newdir12.2)
      ggsave(path = newdir12.2, filename = paste0("PCoA of ", cl, ".2.point.pdf"), p4.2, width = 150, height = 120, units = "mm")
      ggsave(path = newdir12.2, filename = paste0("PCoA of ", cl, ".2.point.png"), p4.2, width = 150, height = 120, units = "mm")
    }
  }
}
