metadata_group <- paste0("metadata_paired.txt")

print("ReporterScore")

library(ggplot2)
library("ReporterScore")

level_s = c("MDD-Pre", "MDD-Post")
pl_list = c("MDD-Pre  MDD-Post")
comparlist = list(c("MDD-Pre", "MDD-Post"))

for (pl in pl_list) {
  pl2 <- as.character(unlist(strsplit(pl, split = "  ")))

  metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  metadata00 <- subset(metadata, Group %in% pl2)
  metadata00$Group <- factor(metadata00$Group,
                             levels = level_s, ordered = TRUE)
  KO_abun <- read.table("out_KO.txt", header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  KO_abun <- KO_abun[, colnames(KO_abun) %in% row.names(metadata00)]
  KO_abun <- KO_abun[which(rowSums(KO_abun) > 0), ]
  metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]

  newdir27 = paste0("ReporterScore/", pl)
  dir.create(newdir27, recursive = TRUE)

  if (length(pl2) == 2) {
    reporter_res <- GRSA(KO_abun, "Group", metadata00,
                         mode = "directed",
                         method = "wilcox.test", perm = 999,
                         type = "pathway", feature = "ko"
    )
  }
  if (length(pl2) != 2) {
    reporter_res <- GRSA(KO_abun, "Group", metadata00,
                         mode = "mixed",
                         method = "kruskal.test", perm = 999,
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

  (p2 <- plot_report_circle_packing(reporter_res, rs_threshold = c(-2.5, 2.5)))
  ggsave(path = newdir27, filename = paste0("report_circle.pdf"), p2, width = 400, height = 350, units = "mm")
  ggsave(path = newdir27, filename = paste0("report_circle.png"), p2, width = 400, height = 350, units = "mm")

  reporter_res_select <- subset(reporter_res$reporter_s, ReporterScore > 2.5 | ReporterScore < -2.5)

  newdir27.1 = paste0("ReporterScore/", pl, "/KOs_in_pathway")
  dir.create(newdir27.1, recursive = TRUE)
  ko_list <- reporter_res_select$ID
  for (ko in ko_list) {
    p3 <- plot_KOs_in_pathway(reporter_res, map_id = ko)
    ggsave(path = newdir27.1, filename = paste0(ko, "_KOs_in_pathway.pdf"), p3, width = 150, height = 150, units = "mm")
    ggsave(path = newdir27.1, filename = paste0(ko, "_KOs_in_pathway.png"), p3, width = 150, height = 150, units = "mm")
  }

  newdir27.2 = paste0("ReporterScore/", pl, "/KOs_distribution")
  dir.create(newdir27.2, recursive = TRUE)
  for (ko in ko_list) {
    p4 <- plot_KOs_distribution(reporter_res, map_id = ko)
    ggsave(path = newdir27.2, filename = paste0(ko, "_KOs_distribution.pdf"), p4, width = 150, height = 150, units = "mm")
    ggsave(path = newdir27.2, filename = paste0(ko, "_KOs_distribution.png"), p4, width = 150, height = 150, units = "mm")
  }

  newdir27.3 = paste0("ReporterScore/", pl, "/KOs_network")
  dir.create(newdir27.3, recursive = TRUE)

  library(myplot)
  for (ko in ko_list) {
    p5 = myplot({
      plot_KOs_network(reporter_res,
                       map_id = ko,
                       main = "", mark_module = TRUE)
    })
    plotsave(p5, path = newdir27.3, file = paste0(ko, "_KOs_network.pdf"), width = 200, height = 200, units = "mm")
    plotsave(p5, path = newdir27.3, file = paste0(ko, "_KOs_network.png"), width = 200, height = 200, units = "mm")
  }

  newdir27.4 = paste0("ReporterScore/", pl, "/KOs_box")
  dir.create(newdir27.4, recursive = TRUE)

  for (ko in ko_list) {
    fit <- try(p6 <- plot_KOs_box(reporter_res, map_id = ko, only_sig = TRUE))
    if ("try-error" %in% class(fit)) {
      next
    } else {
      ggsave(path = newdir27.4, filename = paste0(ko, "_KOs_box.pdf"), p6, width = 350, height = 250, units = "mm")
      ggsave(path = newdir27.4, filename = paste0(ko, "_KOs_box.png"), p6, width = 350, height = 250, units = "mm")
    }
  }

  newdir27.5 = paste0("ReporterScore/", pl, "/KOs_heatmap")
  dir.create(newdir27.5, recursive = TRUE)

  for (ko in ko_list) {
    fit <- try(p7 <- plot_KOs_heatmap(reporter_res, map_id = ko, only_sig = TRUE, heatmap_param = list(cutree_rows = 2)))
    if ("try-error" %in% class(fit)) {
      next
    } else {
      p7 <- plot_KOs_heatmap(reporter_res, map_id = ko, only_sig = TRUE, heatmap_param = list(cutree_rows = 2))
      ggsave(path = newdir27.5, filename = paste0(ko, "_KOs_heatmap.pdf"), p7, width = 350, height = 250, units = "mm")
      ggsave(path = newdir27.5, filename = paste0(ko, "_KOs_heatmap.png"), p7, width = 350, height = 250, units = "mm")
    }
  }

  newdir27.6 = paste0("ReporterScore/", pl, "/KEGG_map")
  dir.create(newdir27.6, recursive = TRUE)

  for (ko in ko_list) {
    fit <- try(plot_KEGG_map(reporter_res, map_id = ko, color_var = "Z_score", save_dir = newdir27.6))
    if ("try-error" %in% class(fit)) {
      next
    } else {
      plot_KEGG_map(reporter_res, map_id = ko, color_var = "Z_score", save_dir = newdir27.6)
    }
  }
}

print("ReporterScore")

library(ggplot2)
library("ReporterScore")

for (pl in pl_list) {
  pl2 <- as.character(unlist(strsplit(pl, split = "  ")))

  metadata = read.table(metadata_group, header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  metadata00 <- subset(metadata, Group %in% pl2)

  KO_abun <- read.table("out_KO.txt", header = T, row.names = 1, sep = "\t", comment.char = "", check.names = F, quot = "", stringsAsFactors = F)
  KO_abun <- KO_abun[, colnames(KO_abun) %in% row.names(metadata00)]
  KO_abun <- KO_abun[which(rowSums(KO_abun) > 0), ]
  metadata00 <- metadata00[order(metadata00$Group, decreasing = TRUE), , drop = FALSE]
  metadata00$Group <- factor(metadata00$Group,
                             levels = level_s, ordered = TRUE)
  newdir27 = paste0("ReporterScore/", pl, "/module/")
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

  (p2 <- plot_report_circle_packing(reporter_res, rs_threshold = c(-2.5, 2.5)))
  ggsave(path = newdir27, filename = paste0("report_circle.pdf"), p2, width = 400, height = 350, units = "mm")
  ggsave(path = newdir27, filename = paste0("report_circle.png"), p2, width = 400, height = 350, units = "mm")

  reporter_res_select <- subset(reporter_res$reporter_s, ReporterScore > 2.5 | ReporterScore < -2.5)

  newdir27.1 = paste0("ReporterScore/", pl, "/module/", "/KOs_in_pathway")
  dir.create(newdir27.1, recursive = TRUE)
  ko_list <- reporter_res_select$ID
  for (ko in ko_list) {
    p3 <- plot_KOs_in_pathway(reporter_res, map_id = ko)
    ggsave(path = newdir27.1, filename = paste0(ko, "_KOs_in_pathway.pdf"), p3, width = 150, height = 150, units = "mm")
    ggsave(path = newdir27.1, filename = paste0(ko, "_KOs_in_pathway.png"), p3, width = 150, height = 150, units = "mm")
  }

  newdir27.2 = paste0("ReporterScore/", pl, "/module/", "/KOs_distribution")
  dir.create(newdir27.2, recursive = TRUE)
  for (ko in ko_list) {
    p4 <- plot_KOs_distribution(reporter_res, map_id = ko)
    ggsave(path = newdir27.2, filename = paste0(ko, "_KOs_distribution.pdf"), p4, width = 150, height = 150, units = "mm")
    ggsave(path = newdir27.2, filename = paste0(ko, "_KOs_distribution.png"), p4, width = 150, height = 150, units = "mm")
  }

  newdir27.3 = paste0("ReporterScore/", pl, "/module/", "/KOs_network")
  dir.create(newdir27.3, recursive = TRUE)

  library(myplot)
  for (ko in ko_list) {
    p5 = myplot({
      plot_KOs_network(reporter_res,
                       map_id = ko,
                       main = "", mark_module = TRUE)
    })
    plotsave(p5, path = newdir27.3, file = paste0(ko, "_KOs_network.pdf"), width = 200, height = 200, units = "mm")
    plotsave(p5, path = newdir27.3, file = paste0(ko, "_KOs_network.png"), width = 200, height = 200, units = "mm")
  }

  newdir27.4 = paste0("ReporterScore/", pl, "/module/", "/KOs_box")
  dir.create(newdir27.4, recursive = TRUE)

  for (ko in ko_list) {
    fit <- try(p6 <- plot_KOs_box(reporter_res, map_id = ko, only_sig = TRUE))
    if ("try-error" %in% class(fit)) {
      next
    } else {
      ggsave(path = newdir27.4, filename = paste0(ko, "_KOs_box.pdf"), p6, width = 350, height = 250, units = "mm")
      ggsave(path = newdir27.4, filename = paste0(ko, "_KOs_box.png"), p6, width = 350, height = 250, units = "mm")
    }
  }

  newdir27.5 = paste0("ReporterScore/", pl, "/module/", "/KOs_heatmap")
  dir.create(newdir27.5, recursive = TRUE)

  for (ko in ko_list) {
    fit <- try(p7 <- plot_KOs_heatmap(reporter_res, map_id = ko, only_sig = TRUE, heatmap_param = list(cutree_rows = 2)))
    if ("try-error" %in% class(fit)) {
      next
    } else {
      p7 <- plot_KOs_heatmap(reporter_res, map_id = ko, only_sig = TRUE, heatmap_param = list(cutree_rows = 2))
      ggsave(path = newdir27.5, filename = paste0(ko, "_KOs_heatmap.pdf"), p7, width = 350, height = 250, units = "mm")
      ggsave(path = newdir27.5, filename = paste0(ko, "_KOs_heatmap.png"), p7, width = 350, height = 250, units = "mm")
    }
  }

  newdir27.6 = paste0("ReporterScore/", pl, "/module/", "/KEGG_map")
  dir.create(newdir27.6, recursive = TRUE)

  for (ko in ko_list) {
    fit <- try(plot_KEGG_map(reporter_res, map_id = ko, color_var = "Z_score", save_dir = newdir27.6))
    if ("try-error" %in% class(fit)) {
      next
    } else {
      plot_KEGG_map(reporter_res, map_id = ko, color_var = "Z_score", save_dir = newdir27.6)
    }
  }
}
