print("ReporterScore")

metadata_group <- "metadata1.txt"
pl_list <- c("HC  MDD")

library(ggplot2)
library("ReporterScore")
set.seed(123)

for (pl in pl_list) {
  pl2 <- as.character(unlist(strsplit(pl, split = "  ")))

  metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  metadata00 <- subset(metadata, Group %in% pl2)

  KO_abun <- read.table("functional_maaslin/out_KO.txt", header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  KO_abun <- KO_abun[, colnames(KO_abun) %in% row.names(metadata00)]
  KO_abun <- KO_abun[which(rowSums(KO_abun) > 0), ]

  newdir27 = paste0("functional_maaslin/reporter_score/", pl, "/module/")
  dir.create(newdir27, recursive = TRUE)

  if (length(pl2) == 2) {
    reporter_res <- GRSA(KO_abun, "Group", metadata00,
                         mode = "directed",
                         method = "wilcox.test", perm = 999,
                         type = "module", feature = "ko"
    )
  }
  if (length(pl2) != 2) {
    reporter_res <- GRSA(KO_abun, "Group", metadata00,
                         mode = "mixed",
                         method = "kruskal.test", perm = 999,
                         type = "module", feature = "ko"
    )
  }
  ko_stat <- reporter_res$ko_stat
  write.table(ko_stat, file = paste0(newdir27, "/ko_stat.txt"), sep = "\t", quote = F, row.names = F)

  reporter_s <- reporter_res$reporter_s
  write.table(reporter_s, file = paste0(newdir27, "/reporter_s.txt"), sep = "\t", quote = F, row.names = F)

  modulelist <- reporter_res$modulelist
  write.table(modulelist, file = paste0(newdir27, "/modulelist.txt"), sep = "\t", quote = F, row.names = F)

  (p1 <- plot_report_bar(reporter_res, rs_threshold = c(-2.5, 2.5), facet_level = TRUE))
  ggsave(path = newdir27, filename = paste0("report_bar.pdf"), p1, width = 370, height = 250, units = "mm")
  ggsave(path = newdir27, filename = paste0("report_bar.png"), p1, width = 370, height = 250, units = "mm")

  (p1 <- plot_report_bar(reporter_res, rs_threshold = c(-2.60, 2.60), facet_level = TRUE))
  ggsave(path = newdir27, filename = paste0("report_bar_filter.pdf"), p1, width = 370, height = 250, units = "mm")
  ggsave(path = newdir27, filename = paste0("report_bar_filter.png"), p1, width = 370, height = 250, units = "mm")

  modulelist <- reporter_res$modulelist
  write.table(modulelist, file = paste0(newdir27, "/modulelist.txt"), sep = "\t", quote = F, row.names = F)
}
