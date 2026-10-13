rm(list = ls())

metadata_group <- "metadata1.txt"
level_s <- c("HC", "MDD")
pl_list <- c("HC  MDD")
comparlist <- list(c("HC", "MDD"))

dir.create("tem", recursive = TRUE)
dir.create("tex", recursive = TRUE)

library(vegan)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(xlsx)

print("alpha diversity index calculation")

newdir3 <- "species_analysis/alpha"
dir.create(newdir3, recursive = TRUE)

for (pl in pl_list) {
  print(pl)

  newdir3 <- paste0("species_analysis/alpha/", pl)
  dir.create(newdir3, recursive = TRUE)

  fit <- try(file.remove(paste0(newdir3, "/alpha.p.xlsx")))
  if ("try-error" %in% class(fit)) {
    next
  } else {
    file.remove(paste0(newdir3, "/alpha.p.xlsx"))
    file.remove(paste0(newdir3, "/alpha.xlsx"))
  }

  cl_list = c("species")
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

pl_list = pl_list

cl_list = c("species")
for (cl in cl_list) {
  print(cl)
  for (pl in pl_list) {
    print(pl)
    pl2 <- as.character(unlist(strsplit(pl, split = "  ")))
    if (length(pl2) != 2) {
      metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
      otutab = read.table(paste0("merged_abundance_table_", cl, ".txt"), header = T, row.names = 1, sep = "\t", quot = "", comment.char = "", check.names = F)
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
      otutab = read.table(paste0("merged_abundance_table_", cl, ".txt"), header = T, row.names = 1, sep = "\t", comment.char = "", quot = "", check.names = F)
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

library("vegan")

metadata <- read.table("metadata_newphenotype2.txt", header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
print("permanova")

continuous <- c("ages", "Waist_Circum", "Hip_Circum", "BMI", "SBP", "DBP", "HAMD", "PHQ", "HAMA", "GAD")
discrete <- c("gender", "Education", "Smoking", "Drinking")

newdir0 = "species_analysis/permanova"
dir.create(newdir0, recursive = TRUE)

if (!is.null(discrete) | !is.null(continuous)) {

  cl_list = c("species")
  kindom_list = c("hebing")
  for (k in kindom_list) {
    for (cl in cl_list) {

      for (pl in pl_list) {

        metadata = read.table("metadata_newphenotype2.txt", header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
        otutab = read.table(paste0("merged_filtered/out_", cl, "_x100_0.01_20filtered.txt"), header = T, row.names = 1, sep = "\t", comment.char = "", quot = "", check.names = F)
        pl2 <- as.character(unlist(strsplit(pl, split = "  ")))
        metadata00 <- subset(metadata, Group %in% pl2)
        metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
        count <- otutab[, colnames(otutab) %in% rownames(metadata00)]
        count <- count[, rownames(metadata00)]

        otu_adonis <- adonis2(
          t(count) ~ .,
          data = metadata00,
          permutations = 9999,
          method = "bray",
          by = "margin"
        )

        print(cl)
        print(otu_adonis$`Pr(>F)`)

        library(ggplot2)

        adonis_df <- data.frame(
          Variable = rownames(otu_adonis),
          R2 = otu_adonis$R2,
          pvalue = otu_adonis$`Pr(>F)`
        )

        adonis_df <- adonis_df[!adonis_df$Variable %in% c("Residual", "Total"), ]

        adonis_df$sig <- cut(
          adonis_df$pvalue,
          breaks = c(-Inf, 0.001, 0.01, 0.05, 1),
          labels = c("***", "**", "*", "")
        )

        adonis_df$color <- ifelse(adonis_df$pvalue < 0.05, "#66619B", "grey70")

        adonis_df$Variable <- factor(adonis_df$Variable, levels = adonis_df$Variable)

        p <- ggplot(adonis_df, aes(x = Variable, y = R2 * 100, fill = color)) +
          geom_col(width = 0.7, show.legend = FALSE) +
          scale_fill_identity() +
          geom_text(aes(label = sig), vjust = 0.20,
                    size = 5) +
          ylab("Explained variance (%)") +
          xlab("") +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.y = element_text(colour = "black"),
            plot.title = element_text(colour = "black")
          )

        print(p)

        ggsave(filename = paste0(newdir0, "/permanova_R2.", cl, "_", k, "_", gsub(" ", "_", pl), ".png"), plot = p, width = 6, height = 4, dpi = 300)
        ggsave(filename = paste0(newdir0, "/permanova_R2.", cl, "_", k, "_", gsub(" ", "_", pl), ".pdf"), plot = p, width = 6, height = 4, dpi = 300)

        library(xlsx)
        write.xlsx(otu_adonis, paste0(newdir0, "/phenotype11.permanova.", cl, "_", k, ".xlsx"),
                   sheetName = pl, append = TRUE,
                   col.names = TRUE,
                   row.names = TRUE)
      }
    }
  }
}

cl_list <- c("bacteria", "fungi", "virus", "archaea")

for (cl in cl_list) {

  data <- read.table(paste0("differential_abundance/maaslin2_filtered/", cl, "/", cl, "_speciesmaaslin_q0.05-5.txt"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")

  combined_matrix <- data
  anno2 <- combined_matrix[, c("feature", "coef")]
  rownames(anno2) <- anno2$feature

  taxnomy <- read.table("taxnomy/taxnomy_zong.txt",
                        sep = "|",
                        na.strings = c("", "NA"),
                        stringsAsFactors = FALSE,
                        header = FALSE)

  taxnomy$V8 <- gsub("s__", "", taxnomy$V8)
  rownames(taxnomy) <- taxnomy$V8
  anno3 <- taxnomy[rownames(anno2), ]

  merge <- merge(anno2, anno3, by = "row.names")

  newdir = paste0(cl, "_itol")
  dir.create(newdir, recursive = TRUE)

  write.table(merge, paste0(newdir, "/MDD_itol_", cl, ".ann.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  color_palette <- c("#BC80BD", "#D9D9D9", "#8DD3C7", "#B2DF8A", "#80B1D3", "#FDB462", "#B3DE69", "#FB8072",
                     "#FCCDE5", "#BEBADA", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#FFFFB3",
                     "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3")

  library(tidyverse)
  color_data <- merge %>%
    mutate(Color = factor(V3, levels = unique(V3), labels = color_palette[1:length(unique(V3))])) %>%
    arrange(V5)

  write.table(color_data, paste0(newdir, "/MDD_itol_", cl, "_maaslin_color.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  write.table(anno3, paste0(newdir, "/MDD_itol_", cl, "_maaslin.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  library(taxonomizr)
  taxa <- as.matrix(read.table(paste0(newdir, "/MDD_itol_", cl, "_maaslin.txt"), sep = "\t"))
  out = makeNewick(taxa)
  cat(out, file = paste0("MDD_itol_", cl, "_maaslin.nwk"))

}

cl_list <- c("bacteria", "fungi", "virus", "archaea")

for (cl in cl_list) {

  B <- read.table(paste0("differential_abundance/maaslin2_filtered/", "bacteria", "/", "bacteria", "_speciesmaaslin_q0.05-5.txt"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")
  F <- read.table(paste0("differential_abundance/maaslin2_filtered/", "fungi", "/", "fungi", "_speciesmaaslin_q0.05-5.txt"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")
  V <- read.table(paste0("differential_abundance/maaslin2_filtered/", "virus", "/", "virus", "_speciesmaaslin_q0.05-5.txt"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")
  A <- read.table(paste0("differential_abundance/maaslin2_filtered/", "archaea", "/", "archaea", "_speciesmaaslin_q0.05-5.txt"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")

  B$Kingdom <- "bacteria"
  F$Kingdom <- "fungi"
  V$Kingdom <- "virus"
  A$Kingdom <- "archaea"

  data <- rbind(B, F, V, A)

  data$Taxon <- rownames(data)

  combined_matrix <- data
  anno2 <- combined_matrix[, c("feature", "coef")]
  rownames(anno2) <- anno2$feature

  taxnomy <- read.table("taxnomy/taxnomy_zong.txt",
                        sep = "|",
                        na.strings = c("", "NA"),
                        stringsAsFactors = FALSE,
                        header = FALSE)

  taxnomy$V8 <- gsub("s__", "", taxnomy$V8)
  rownames(taxnomy) <- taxnomy$V8
  anno3 <- taxnomy[rownames(anno2), ]

  combined <- rbind(B, F, V, A)

  anno_tax <- taxnomy[combined$feature, ]

  final_merge <- cbind(combined, anno_tax)

  str(final_merge)
  head(final_merge)

  merge <- merge(anno2, anno3, by = "row.names")
  cl = "hebing"
  newdir = paste0(cl, "_itol")
  dir.create(newdir, recursive = TRUE)
  write.table(final_merge, paste0(newdir, "/MDD_itol_hebing_maaslin_tax.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  write.table(merge, paste0(newdir, "/MDD_itol_", cl, ".ann.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  color_palette <- c("#BC80BD", "#D9D9D9", "#8DD3C7", "#B2DF8A", "#80B1D3", "#FDB462", "#B3DE69", "#FB8072",
                     "#FCCDE5", "#BEBADA", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#FFFFB3",
                     "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3")

  library(tidyverse)
  color_data <- merge %>%
    mutate(Color = factor(V3, levels = unique(V3), labels = color_palette[1:length(unique(V3))])) %>%
    arrange(V5)

  write.table(color_data, paste0(newdir, "/MDD_itol_", cl, "_maaslin_color.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  write.table(anno3, paste0(newdir, "/MDD_itol_", cl, "_maaslin.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  library(taxonomizr)
  taxa <- as.matrix(read.table(paste0(newdir, "/MDD_itol_", cl, "_maaslin.txt"), sep = "\t"))
  out = makeNewick(taxa)
  cat(out, file = paste0("MDD_itol_", cl, "_maaslin.nwk"))

}

library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(here)
library(PupillometryR)
library(readxl)

species <- read_xlsx("selected_species.xlsx", sheet = "Sheet1")

species_top20 <- species$feature

otutab <- read.table(
  "prediction_maaslin/merge_mutil_kindom_change.txt",
  header = TRUE, row.names = 1, sep = "\t",
  stringsAsFactors = FALSE, check.names = FALSE
)

metadata <- read.table(
  "metadata1.txt",
  header = TRUE, row.names = 1, sep = "\t",
  stringsAsFactors = FALSE
)

otutab_top20 <- otutab[species_top20, ]

plot_matrix <- as.data.frame(t(otutab_top20))

plot_matrix$Group <- metadata[rownames(plot_matrix), "Group"]

table(plot_matrix$Group)

level_s <- c("HC", "MDD")
comparlist <- list(c("HC", "MDD"))

plot_matrix$Group <- factor(
  plot_matrix$Group,
  levels = level_s,
  ordered = TRUE
)
index_list <- colnames(plot_matrix)[!colnames(plot_matrix) == "Group"]

outdir <- "boxplot"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (index1 in index_list) {

  df <- plot_matrix[, c(index1, "Group"), drop = FALSE]
  colnames(df) <- c("value", "Group")

  p2 <- ggplot(df, aes(Group, log(value), fill = Group)) +
    geom_flat_violin(aes(fill = Group), position = position_nudge(x = .25), color = "black", alpha = 0.7) +
    geom_jitter(position = position_jitter(0.15), size = 2.5, alpha = 0.7, aes(color = Group), na.rm = TRUE, inherit.aes = TRUE) +
    geom_boxplot(width = .1, position = position_nudge(x = 0.25), fill = "white", size = 0.5) +

    labs(x = "", y = paste0(index1)) +
    geom_signif(
      comparisons = comparlist,
      na.rm = TRUE,
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05, " " = 2),
      test = wilcox.test,
      step_increase = 0.1,
      tip_length = 0,
      size = 0.9,
      col = "black",
      textsize = 10,
      vjust = 0.68,
      position = "identity"

    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12, face = "bold", colour = 'black'),
      plot.title = element_text(size = 15, face = "bold", vjust = 0.5, hjust = 0.5, colour = 'black'),
      legend.position = "bo",
      text = element_text(family = "sans", size = 10, colour = 'black'),
      axis.text.y = element_text(family = "sans", size = 8, colour = 'black'),
      axis.text.x = element_text(size = 12, hjust = 0.5, face = "bold", colour = "black",
                                 vjust = 0.5, family = "sans"))

  p2
  p2 = p2 +
    scale_color_manual(values = c("HC" = "#6883C8", "MDD" = "#EB6468")) +
    scale_fill_manual(values = c("HC" = "#6883C8", "MDD" = "#EB6468"))
  newdir = "boxplot"
  dir.create(newdir, recursive = TRUE)
  ggsave(path = newdir, filename = paste0(index1, ".2.pdf"), p2, width = 80, height = 95, units = "mm")
  ggsave(path = newdir, filename = paste0(index1, ".2.png"), p2, width = 80, height = 95, units = "mm")

}

library(tidyverse)
library(ggcor)
library(vegan)
library(readxl)
library(openxlsx)

pheno <- read_xlsx("phenotype_mantel/phenotype_data.xlsx", sheet = "Sheet3") %>%
  column_to_rownames("ID")

pheno_num <- pheno %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.))))

microbe <- read.table(
  "merged_filtered/out_species_x100_0.01_20filtered.txt",
  header = TRUE, sep = "\t", check.names = FALSE
) %>%
  column_to_rownames("ID")

tax <- read.table(
  "taxnomy/taxnomy_zong.txt",
  sep = "|",
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE,
  header = FALSE,
  comment.char = ""
) %>%
  setNames(c("Kingdom", "Kingdom2", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

tax <- tax %>%
  mutate(
    ID = gsub("^s__", "", Species),
    Kingdom = gsub("^k__", "", Kingdom)
  ) %>%
  select(ID, Kingdom)

rownames(tax) <- tax$ID

group <- read.table("metadata1.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

common_samples <- intersect(colnames(microbe), rownames(pheno_num))

microbe <- microbe[, common_samples, drop = FALSE]
pheno_num <- pheno_num[common_samples, , drop = FALSE]

group <- group %>%
  filter(ID %in% common_samples) %>%
  arrange(match(ID, common_samples))

tax <- tax %>%
  filter(ID %in% rownames(microbe)) %>%
  arrange(match(ID, rownames(microbe)))

microbe <- microbe[rownames(microbe) %in% tax$ID, , drop = FALSE]

microbe_t <- t(microbe) %>% as.data.frame()
microbe_num <- microbe_t %>% mutate(across(everything(), ~ as.numeric(.)))

spec_list <- tax %>%
  select(ID, Kingdom) %>%
  mutate(Kingdom = coalesce(Kingdom, "Unclassified")) %>%
  mutate(Kingdom = gsub("^d__", "", Kingdom)) %>%
  group_by(Kingdom) %>%
  summarise(
    otus = list(match(ID, rownames(microbe)))
  ) %>%
  deframe() %>%
  compact() %>%
  keep(~ length(.x) >= 5)

mantel_raw <- mantel_test(
  microbe_num, pheno_num,
  spec.select = spec_list,
  mantel.fun = "mantel",
  method = "spearman",
  permutations = 999

)

mantel <- mantel_raw %>%
  mutate(
    rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
             labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
             labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
    pcol = case_when(
      p.value >= 0.05 ~ "p >= 0.05",
      p.value < 0.05 & r > 0 ~ "p < 0.05, r > 0",
      p.value < 0.05 & r < 0 ~ "p < 0.05, r < 0",
      TRUE ~ NA_character_
    )
  )

mantel_output <- mantel_raw %>%
  select(spec = spec, env = env, mantel_r = r, p_value = p.value)

write.xlsx(mantel_output, "mantel_test_results21.xlsx", rowNames = FALSE)

quickcor(pheno_num, type = "upper", show.diag = FALSE, side = "top") +
  geom_circle2() +

  anno_link(
    aes(
      colour = pcol,
      size = abs(r)
    ),
    data = mantel,
    curvature = 0.2
  ) +

  scale_fill_gradient2(
    low = "#377EB8",
    mid = "#F7F7F7",
    high = "#E41A1C",
    midpoint = 0
  ) +

  scale_colour_manual(
    name = "Mantel test",
    values = c(
      "p >= 0.05"        = "grey70",
      "p < 0.05, r > 0" = "#E41A1C",
      "p < 0.05, r < 0" = "#377EB8"
    ),
    labels = c(
      "p >= 0.05"        = "p >= 0.05 (ns)",
      "p < 0.05, r > 0" = "p < 0.05, r > 0",
      "p < 0.05, r < 0" = "p < 0.05, r < 0"
    )
  ) +

  scale_size_continuous(
    name = "|Mantel's r|",
    range = c(0.1, 1),
    limits = c(0, 0.2)
  ) +

  guides(
    colour = guide_legend(title = "Mantel's p", order = 1),
    size   = guide_legend(title = "|Mantel's r|", order = 2),
    fill   = guide_colorbar(title = "Spearman's r", order = 3)
  ) +

  theme(
    axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 9),
    axis.text.y.right = element_text(size = 9),
    axis.text.y.left = element_blank(),
    axis.ticks.y.left = element_blank()
  )

ggsave("style1_minimal_circle21.png", height = 10, width = 14)
ggsave("style1_minimal_circle21.pdf", height = 10, width = 14)
