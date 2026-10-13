rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tidyverse)
library(pheatmap)

feature_outdir = "anxiety_code_Results/"

cl_list <- c("B", "F", "M", "F_B", "B_M", "F_M", "B_F_M")
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

hebing$Group <- factor(hebing$Group, levels = c("B", "F", "M", "F_B", "B_M", "F_M", "B_F_M"))

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

cl_list <- c("B", "F", "M", "F_B", "B_M", "F_M", "B_F_M")

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

hebing$Group <- factor(hebing$Group,
                       levels = c("B", "F", "M", "F_B", "B_M", "F_M", "B_F_M"), ordered = TRUE)

first_row_df <- hebing %>%
  group_by(Group) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))

first_row_df$Group <- factor(first_row_df$Group,
                             levels = c("B", "F", "M", "F_B", "B_M", "F_M", "B_F_M"), ordered = TRUE)

write.csv(first_row_df, file = "anxiety_Discovery_selfcv_repeat_max_ALL_by_group.csv", quote = F, row.names = F)

p_Sensitivity <- ggplot(data = first_row_df, aes(x = Group, y = Sensitivity)) +
  geom_col(fill = "#1F78B4", color = "black") +
  geom_text(aes(label = percent(Sensitivity, accuracy = 0.1)), size = 6, vjust = 0.5, color = "white", hjust = 1) +
  theme_classic() +
  ggtitle("Sensitivity") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"))

p_Sensitivity

p_Specificity <- ggplot(data = first_row_df, aes(x = Group, y = Specificity)) +
  geom_col(fill = "#33A02C", color = "black") +
  geom_text(aes(label = percent(Specificity, accuracy = 0.1)), vjust = 0.5, size = 6, color = "white", hjust = 1) +
  theme_classic() +
  ggtitle("Specificity") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"))

p_Specificity

p_AUC <- ggplot(data = first_row_df, aes(x = Group, y = AUC)) +
  geom_col(fill = "#BC80BD", color = "black") +
  geom_text(aes(label = percent(AUC, accuracy = 0.1)), vjust = 0.5, size = 6, color = "white", hjust = 1) +
  theme_classic() +
  ggtitle("AUC") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"))

p_AUC

p_Accuracy <- ggplot(data = first_row_df, aes(x = Group, y = Balanced.Accuracy)) +
  geom_col(fill = "grey50", color = "black") +
  geom_text(aes(label = percent(Balanced.Accuracy, accuracy = 0.1)), vjust = 0.5, size = 6, color = "white", hjust = 1) +
  theme_classic() +
  ggtitle("Accuracy") +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, color = "black"))

p_Accuracy

color_palette <- c("#BC80BD", "#D9D9D9", "#8DD3C7", "#B2DF8A", "#80B1D3", "#FDB462", "#B3DE69", "#FB8072",
                   "#FCCDE5", "#BEBADA", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#FFFFB3",
                   "#E31A1C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3")

plot_df <- first_row_df %>%
  as.data.frame() %>%
  transmute(
    Group = factor(Group, levels = rev(c("B", "F", "M", "F_B", "B_M", "F_M", "B_F_M")), ordered = TRUE),
    Sensitivity = Sensitivity,
    Specificity = Specificity,
    Accuracy    = Balanced.Accuracy
  ) %>%
  pivot_longer(-Group, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Sensitivity", "Specificity", "Accuracy")))

p <- ggplot(plot_df, aes(x = Group, y = Value)) +
  geom_col(aes(fill = Metric), color = "black", width = 0.7, show.legend = FALSE) +
  geom_text(aes(label = percent(Value, accuracy = 0.1)),
            size = 4.5, color = "white", hjust = 1.2) +
  scale_fill_manual(values = c(Sensitivity = "#80B1D3",
                               Specificity = "#5FB399",
                               Accuracy    = "#BC80BD")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_flip() +
  facet_grid(. ~ Metric) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(size = 15, color = "black"),
    strip.background = element_blank(),
    axis.text.y = element_text(size = 13, color = "black"),
    axis.text.x = element_blank(),
    axis.line.x = element_blank(), axis.ticks.x = element_blank(),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

ggsave("combined_plot2.pdf", p, width = 130, height = 120, units = "mm")
ggsave("combined_plot2.png", p, width = 130, height = 120, units = "mm")

df_wide_new1 <- as.matrix(read.table(paste0("feature_all.txt"), header = T, row.names = 1, sep = '\t', quot = "", check.names = F))

colnames <- colnames(df_wide_new1)

convert_to_group <- function(colname) {
  num_parts <- length(strsplit(colname, "_")[[1]])

  if (num_parts == 1) {
    return("one")
  } else if (num_parts == 2) {
    return("two")
  } else if (num_parts == 3) {
    return("three")
  } else {
    return(NA_character_)
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
M <- read.table(paste0("M_significant_results.tsv"), header = T, row.names = 1, sep = '\t', comment.char = "", quot = "")

combined_matrix <- rbind(B, F, M)
combined_matrix$feature <- rownames(combined_matrix)
anno2 <- combined_matrix[rownames(df_wide_new1), ][, c("feature", "coef")]

taxnomy <- read.table(paste0("../taxnomy/taxnomy_zong.txt"),
                      sep = "|",
                      na.strings = c("", "NA"),
                      stringsAsFactors = FALSE,
                      header = FALSE)

taxnomy$feature <- gsub("s__", "", taxnomy$V8)
rownames(taxnomy) <- taxnomy$feature

features <- rownames(df_wide_new1)

idx <- match(features, rownames(taxnomy))

anno_V1 <- rep("Amino_acid", length(features))

anno_V1[!is.na(idx)] <- taxnomy$V1[idx[!is.na(idx)]]

anno3 <- data.frame(
  feature = features,
  V1      = anno_V1,
  stringsAsFactors = FALSE
)

rownames(anno3) <- features

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
                                 levels = c("Bacteria", "Fungi", "Amino_acid"), ordered = TRUE)

colnames(annotation_row1)

ann_colors = list(
  Group = c(one = "#3E4989", two = "#A6CEE3", three = "#CCEBC5"),
  Kindom = c(Bacteria = "#7570B3", Fungi = "#BC80BD", Amino_acid = "#FDDAEC"),
  Change = c(Up = "#DE77AE", Down = "#7FBC41"),
  importance = c(colorRampPalette(c('#3B4992FF', 'white'))(60),
                 colorRampPalette(c('white', 'white'))(2),
                 colorRampPalette(c('white', '#EE0000FF'))(60)))
mycol <- c(colorRampPalette(c("#A6CEE3", "#3B4992FF"))(100))

df_wide_new2 <- df_wide_new1[rownames(annotation_row2), ]
df_wide_new2[is.na(df_wide_new2)] <- 0

p3 = pheatmap(t(df_wide_new2),
              cellwidth = 6, cellheight = 16,
              cluster_cols = FALSE,
              cluster_rows = FALSE,
              color = mycol,
              fontsize_number = 12,
              number_color = "white",
              fontsize_row = 11,
              fontsize_col = 7,
              annotation_row = annotation_row1,
              annotation_col = annotation_row2,
              annotation_colors = ann_colors
)

ggsave(filename = paste0("heatmap.feature.pdf"), p3, width = 300, height = 150, units = "mm", limitsize = FALSE)
ggsave(filename = paste0("heatmap.feature.png"), p3, width = 300, height = 150, units = "mm", limitsize = FALSE)
