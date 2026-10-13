library(ggplot2)
library(tidyverse)
library(pheatmap)

feature_outdir = "anxiety_code_Results/"

cl_list <- c("B", "F", "V", "A",
             "A_B", "B_V", "A_V", "F_B", "F_V", "A_F",
             "A_B_F", "A_B_V", "A_F_V", "F_B_V",
             "A_F_B_V")

hebing <- NULL

for (cl in cl_list) {
  file <- read.csv(paste0(feature_outdir, cl, "/anxiety_Discovery_selfcv_repeat_max_", cl, ".csv"), row.names = 1)
  file$Group <- cl
  if (is.null(hebing)) {
    hebing <- file
  } else {
    hebing <- rbind(hebing, file)
  }
}

write.csv(hebing, file = "anxiety_Discovery_selfcv_repeat_max_ALL.csv", quote = F, row.names = F)

hebing$Group <- factor(hebing$Group, levels = c("B", "F", "V", "A",
                                                "A_B", "B_V", "A_V", "F_B", "F_V", "A_F",
                                                "A_B_F", "A_B_V", "A_F_V", "F_B_V",
                                                "A_F_B_V", "Top30"))

col = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
        "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#FF7F00", "#EF5656")

p1 <- ggplot(hebing, aes(x = Group, y = AUC, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.15) +
  geom_boxplot(size = 0.5, fill = "white", outlier.shape = NA) +
  geom_jitter(aes(fill = Group), width = 0.2, shape = 21, size = 2.5) +
  theme_bw() +
  scale_fill_manual(values = col) +
  ggtitle("") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 13, color = "black"),
        axis.title.y = element_text(size = 15, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = FALSE)
p1

ggsave(filename = paste0("boxplot_top.pdf"), p1, width = 150, height = 110, units = "mm")
ggsave(filename = paste0("boxplot_top.png"), p1, width = 150, height = 110, units = "mm")

df_wide_new1 <- as.matrix(read.table(paste0("feature_all.txt"), header = T, row.names = 1, sep = '\t', quot = "", check.names = F))

colnames <- colnames(df_wide_new1)

convert_to_group <- function(colname) {
  num_letters <- nchar(colname)

  if (num_letters == 1) {
    return("one")
  } else if (num_letters == 2) {
    return("two")
  } else if (num_letters == 3) {
    return("three")
  } else if (num_letters == 4) {
    return("four")
  }
}

Group <- sapply(colnames, convert_to_group)

annotation_row1 <- data.frame(ID = colnames, Group = Group)
annotation_row1 <- annotation_row1 %>% select(-ID)
annotation_row1$Group <- factor(annotation_row1$Group,
                                levels = c("one", "two", "three", "four"), ordered = TRUE)

df_wide_new1[is.na(df_wide_new1)] <- 0
row_sums <- rowSums(df_wide_new1)

anno1 <- data.frame(
  feature = rownames(df_wide_new1),
  RowSum = row_sums
)

B <- read.table(paste0("B_significant_results.tsv"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")
F <- read.table(paste0("F_significant_results.tsv"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")
V <- read.table(paste0("V_significant_results.tsv"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")
A <- read.table(paste0("A_significant_results.tsv"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")

combined_matrix <- rbind(B, F, V, A)
combined_matrix$feature <- rownames(combined_matrix)
anno2 <- combined_matrix[rownames(df_wide_new1), ][, c("feature", "coef")]

taxnomy <- read.table(paste0("../taxnomy/taxnomy_zong.txt"),
                      sep = "|",
                      na.strings = c("", "NA"),
                      stringsAsFactors = FALSE,
                      header = FALSE)

taxnomy$feature <- gsub("s__", "", taxnomy$V8)
rownames(taxnomy) <- taxnomy$feature
anno3 <- taxnomy[rownames(df_wide_new1), c("feature", "V1")]

merged_matrix <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), list(anno1, anno2, anno3))

print(merged_matrix)

rownames(merged_matrix) <- merged_matrix$feature

annotation_row2 <- merged_matrix %>%
  mutate(
    Change = ifelse(coef > 0, "Up", "Down"),
    V1 = sub("^d__", "", V1),
    V1 = gsub("Eukaryota", "Fungi", V1),
    importance = rank(RowSum, ties.method = "first")
  ) %>%
  rename(
    Kindom = V1
  )
write.table(annotation_row2, paste0("anno.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)

annotation_row2 <- annotation_row2 %>%
  select(-coef, -feature, -RowSum) %>%
  arrange(importance)

annotation_row2$Kindom <- factor(annotation_row2$Kindom,
                                 levels = c("Bacteria", "Viruses", "Archaea", "Fungi"), ordered = TRUE)

colnames(annotation_row1)

ann_colors = list(
  Group = c(one = "#3E4989", two = "#386CB0", three = "#A6CEE3", four = "#CCEBC5"),
  Kindom = c(Bacteria = "#7570B3", Viruses = "#BC80BD", Archaea = "#DECBE4", Fungi = "#FDDAEC"),
  Change = c(Up = "#DE77AE", Down = "#7FBC41"),
  importance = c(colorRampPalette(c('#3B4992FF', 'white'))(60),
                 colorRampPalette(c('white', 'white'))(2),
                 colorRampPalette(c('white', '#EE0000FF'))(60)))
mycol <- c(colorRampPalette(c("#A6CEE3", "#3B4992FF"))(100))

df_wide_new2 <- df_wide_new1[rownames(annotation_row2), ]
df_wide_new2[is.na(df_wide_new2)] <- 0

p3 = pheatmap(t(df_wide_new2),
              cellwidth = 8, cellheight = 10,
              cluster_cols = FALSE,
              cluster_rows = FALSE,
              color = mycol,
              fontsize_number = 12,
              number_color = "white",
              fontsize_row = 11,
              fontsize_col = 8,
              annotation_row = annotation_row1,
              annotation_col = annotation_row2,
              annotation_colors = ann_colors
)

ggsave(filename = paste0("heatmap.feature.pdf"), p3, width = 300, height = 150, units = "mm", limitsize = FALSE)
ggsave(filename = paste0("heatmap.feature.png"), p3, width = 300, height = 150, units = "mm", limitsize = FALSE)
