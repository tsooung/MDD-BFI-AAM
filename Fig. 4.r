print("Fig. 4A OPLS-DA scores plot of serum metabolic profiles, MDD versus HC")

suppressMessages(library(tidyverse))
suppressMessages(library(MetaboAnalystR))
library(openxlsx)
library(RColorBrewer)
library(ggrepel)

infile <- "metabolite_matrix.tsv"
group_file <- "sample_group.txt"
outdir <- "."
prefix <- "result"
width <- 9.6
height <- 8

outPrefix <- paste0(outdir, "/", prefix)

inDF <- read_tsv(infile, col_names = T,
                 col_types = cols(), locale = locale(encoding = "UTF-8"))
gDF <- read_tsv(group_file, col_types = "cc",
                locale = locale(encoding = "UTF-8"), na = "") %>%
  mutate(Group = factor(Group, levels = unique(Group)))

groupNum <- length(unique(gDF$Group))
colorSets <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#377eb8", "#e41a1c", "#984ea3", "#666666")
colorSets <- colorSets[1:groupNum]
names(colorSets) <- levels(gDF$Group)
colorSets1 <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), brewer.pal(12, "Set3"))
shapeSets <- c(rep(21, groupNum))

inDFReplace <- inDF %>%
  select(c("Index", gDF$Sample)) %>%
  mutate(IndexReplace = str_c("Metabo", 1:nrow(inDF))) %>%
  select(IndexReplace, everything())
inDF <- inDFReplace %>% select(-Index)
cDF <- inDF %>%
  gather(key = Sample, value = value, 2:ncol(inDF)) %>%
  spread(key = names(inDF)[1], value = "value")
matDF <- gDF %>% inner_join(cDF, by = "Sample")

opls_log <- "LogNorm"
opls_scale <- "MeanCenter"

mSet <- InitDataObjects("conc", "stat", FALSE)
mSet <- Read.TextData(mSet, "data4PLS.csv", "rowu", "disc")
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)
mSet <- Normalization(mSet, "NULL", opls_log, opls_scale, "S10T0", ratio = FALSE, ratioNum = 20)
mSet <- OPLSR.Anal(mSet, reg = TRUE)
orthoScore <- mSet$analSet$oplsda$orthoScoreMN
score <- mSet$analSet$oplsda$scoreMN
scoreDF <- score %>% cbind(orthoScore) %>%
  as.data.frame %>% rownames_to_column("Sample") %>%
  as_tibble %>% inner_join(gDF, by = "Sample")
p1 <- mSet$analSet$oplsda$modelDF$R2X[1]
o1 <- mSet$analSet$oplsda$modelDF$R2X[2]

g <- ggplot(scoreDF, aes(x = p1, y = o1, colour = Group, label = Sample)) +
  geom_point(aes(fill = Group, shape = Group, colour = Group)) +
  scale_color_manual(values = colorSets) +
  scale_fill_manual(values = colorSets) +
  scale_shape_manual(values = shapeSets) +
  theme_bw() +
  labs(title = "Scores OPLS-DA Plot",
       x = paste("T score[1] ( ", p1 * 100, "% )", sep = ""),
       y = paste("Orthogonal T score[1] ( ", o1 * 100, "% )", sep = "")) +
  theme(plot.title = element_text(hjust = .5))
gHtml <- g
g <- g + theme(legend.position = "right",
               legend.box = "horizontal",
               legend.title = element_text(face = "bold"))

g <- g + stat_ellipse(show.legend = F, geom = "polygon", alpha = 0.2, aes(fill = Group))
gHtml <- gHtml + stat_ellipse(show.legend = F, geom = "polygon", alpha = 0.2, aes(fill = Group))

if (nrow(gDF) <= 20) {
  set.seed(100)
  g <- g + geom_label_repel(show.legend = F, size = 2.5, segment.size = .4,
                            segment.alpha = .4, box.padding = .6, point.padding = .5,
                            arrow = arrow(length = unit(0.01, "npc"),
                                          type = "closed", ends = "first"))
}

png(file = paste(outPrefix, "_OPLS-DA_scorePlot.png", sep = ""),
    width = width, height = height, res = 300, units = "in")
print(g)
dev.off()
pdf(file = paste(outPrefix, "_OPLS-DA_scorePlot.pdf", sep = ""),
    width = width, height = height)
print(g)
dev.off()

mSet <- PlotOPLS.Splot(mSet, "sPlot", format = "png", dpi = 300, width = NA, "all")
oplsVIP <- mSet$analSet$oplsda$vipVn %>%
  as.data.frame %>%
  `colnames<-`("VIP") %>%
  rownames_to_column("IndexReplace") %>%
  as_tibble %>%
  inner_join(inDFReplace[c("IndexReplace", "Index")], by = "IndexReplace") %>%
  select(Index, everything())

sPlotDF <- mSet$analSet$oplsda$splot.mat %>%
  as.data.frame %>%
  rownames_to_column("IndexReplace") %>%
  as_tibble %>%
  inner_join(oplsVIP, by = "IndexReplace") %>%
  mutate(Color = case_when(VIP > 1 ~ "#f80903", TRUE ~ "#12bd0f"))

oplsVIP <- oplsVIP %>% select(-IndexReplace)
write_tsv(oplsVIP, paste0(outPrefix, "_OPLS-DA_VIP.txt"))

g <- ggplot(sPlotDF, aes(x = `p[1]`, y = `p(corr)[1]`, label = Index, text = VIP)) +
  geom_hline(yintercept = 0, colour = "gray", size = 0.75) +
  geom_vline(xintercept = 0, colour = "gray", size = 0.75) +
  geom_jitter(colour = sPlotDF$Color, alpha = 0.7) +
  theme_bw() +
  labs(title = "OPLS-DA S-Plot") +
  theme(plot.title = element_text(hjust = .5))
png(file = paste(outPrefix, "_OPLS-DA_SPlot.png", sep = ""),
    width = width, height = height, res = 300, units = "in")
print(g)
dev.off()
pdf(file = paste(outPrefix, "_OPLS-DA_SPlot.pdf", sep = ""),
    width = width, height = height)
print(g)
dev.off()

modelDF <- mSet$analSet$oplsda$modelDF %>%
  rownames_to_column("Type") %>%
  select(one_of("Type", "R2X", "R2Y", "Q2")) %>%
  slice(1:2) %>% gather(key = Param, value = Value, R2X, R2Y, Q2) %>%
  mutate(Param = factor(Param, levels = c("R2X", "R2Y", "Q2")))

g <- ggplot(modelDF, aes(x = Type, y = Value, fill = Param, label = Value)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(vjust = -0.15, position = position_dodge(0.9)) +
  scale_fill_manual(name = "",
                    values = colorSets1,
                    labels = c(expression(paste("R"^{2}, "X")),
                               expression(paste("R"^{2}, "Y")),
                               expression(paste("Q"^{2})))) +
  scale_x_discrete(limits = c("p1", "o1")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(legend.text = element_text(hjust = 0))
g <- g + theme(legend.position = c(1, 1),
               legend.justification = c(1, 1),
               legend.background = element_rect(fill = alpha("white", 0.3),
                                                colour = "black", size = 0.3))
png(file = paste(outPrefix, "_OPLS-DA_model.png", sep = ""),
    width = width, height = height, res = 300, units = "in")
print(g)
dev.off()
pdf(file = paste(outPrefix, "_OPLS-DA_model.pdf", sep = ""),
    width = width, height = height)
print(g)
dev.off()

PlotOPLS.Permutation <- function(mSetObj = NA, imgName, num = 100, width = NA, height = NA) {
  mSetObj <- .get.mSet(mSetObj)
  if (mSetObj$analSet$oplsda$reg) {
    cls <- scale(as.numeric(mSetObj$dataSet$cls))[, 1]
  } else {
    cls <- model.matrix(~ mSetObj$dataSet$cls - 1)
  }

  datmat <- as.matrix(mSetObj$dataSet$norm)

  cv.num <- min(7, dim(mSetObj$dataSet$norm)[1] - 1)
  perm.res <- perform_opls(datmat, cls, predI = 1, orthoI = NA, permI = num, crossvalI = cv.num)
  r.vec <- perm.res$suppLs[["permMN"]][, "R2Y(cum)"]
  q.vec <- perm.res$suppLs[["permMN"]][, "Q2(cum)"]

  better.rhits <- sum(r.vec[-1] > r.vec[1])

  if (better.rhits == 0) {
    pr <- paste("p < ", 1 / num, " (", better.rhits, "/", num, ")", sep = "")
  } else {
    p <- better.rhits / num
    pr <- paste("p = ", signif(p, digits = 5), " (", better.rhits, "/", num, ")", sep = "")
  }
  better.qhits <- sum(q.vec[-1] > q.vec[1])
  if (better.qhits == 0) {
    pq <- paste("p < ", 1 / num, " (", better.qhits, "/", num, ")", sep = "")
  } else {
    p <- better.qhits / num
    pq <- paste("p = ", signif(p, digits = 5), " (", better.qhits, "/", num, ")", sep = "")
  }

  rng <- range(c(r.vec, q.vec, 1))

  suppressMessages(library(tidyverse))
  library(RColorBrewer)
  colorSets <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Accent"), brewer.pal(12, "Set3"))
  r2y = rep("Perm R2Y", 200)
  q2 = rep("Perm Q2", 200)

  rawQ2 <- mSetObj$analSet$oplsda$summaryDF[["Q2(cum)"]]
  rawR2Y <- mSetObj$analSet$oplsda$summaryDF[["R2Y(cum)"]]
  rawR2X <- mSetObj$analSet$oplsda$summaryDF[["R2X(cum)"]]

  oplsPermutation <- cbind(basename(getwd()), rawR2Y, rawQ2, better.rhits / num, better.qhits / num)
  colnames(oplsPermutation) <- c("groupName", "R2", "Q2", "pr", "pq")
  oplsPermutation <- as.data.frame(oplsPermutation)
  write.csv(oplsPermutation, row.names = F, file = "oplsda_permutation.csv")

  permDF = tibble(`Values` = c(r.vec[-1], q.vec[-1]), `model_type` = c(r2y, q2))
  g <- ggplot(permDF, aes(x = Values, fill = model_type)) +
    geom_histogram(bins = 30, color = "white", size = 0.2) +
    scale_fill_manual(values = c(colorSets[3], colorSets[2]),
                      labels = c(expression(paste("Perm ", "Q"^{2}, "  ")), expression(paste("Perm ", "R"^{2}, "Y")))) +
    theme(legend.text = element_text(hjust = 0))
  h <- max(ggplot_build(g)$data[[1]][["count"]])
  g <- g + annotate("segment", x = rawQ2, xend = rawQ2, y = h / 2, yend = 0, size = 0.7, arrow = arrow(length = unit(0.01, "npc"))) +
    annotate("text", x = rawQ2, y = h / 1.8, label = c(substitute(Q^2 ~ "=" ~ a, list(a = rawQ2))), parse = TRUE, size = 3, hjust = 0.8, fontface = "plain") +
    annotate("text", x = rawQ2, y = h / 1.9, label = pq, size = 3, hjust = 0.8, fontface = "plain") +
    annotate("segment", x = rawR2Y, xend = rawR2Y, y = h / 3, yend = 0, size = 0.7, arrow = arrow(length = unit(0.01, "npc"))) +
    annotate("text", x = rawR2Y, y = h / 2.5, label = c(substitute(R^2 ~ "Y" ~ "=" ~ a, list(a = rawR2Y))), parse = TRUE, size = 3, hjust = 0.8, fontface = "plain") +
    annotate("text", x = rawR2Y, y = h / 2.7, label = pr, size = 3, hjust = 0.8, fontface = "plain") +
    annotate("segment", x = rawR2X, xend = rawR2X, y = h / 4, yend = 0, size = 0.7, arrow = arrow(length = unit(0.01, "npc"))) +
    annotate("text", x = rawR2X, y = h / 3.5, label = c(substitute(R^2 ~ "X" ~ "=" ~ a, list(a = rawR2X))), parse = TRUE, size = 3, hjust = 0.8, fontface = "plain") +
    theme_classic() + labs(x = "Permutations", y = "Frequency") +
    scale_x_continuous(breaks = seq(-3, 1, 0.2)) +
    theme(legend.title = element_blank(), legend.position = c(0, 1), legend.justification = c(0, 1),
          legend.background = element_rect(fill = alpha("white", 0), colour = "black", size = 0.15),
          panel.grid.minor.x = element_blank())

  png(file = paste(imgName, ".png", sep = ""), width = width, height = height, res = 300, units = "in")
  print(g)
  dev.off()
  pdf(file = paste(imgName, ".pdf", sep = ""), width = width, height = height)
  print(g)
  dev.off()
  mSetObj$imgSet$opls.permut <- imgName

  msg <- paste("Empirical p-values Q2: ", pq, " and R2Y: ", pr)
  if (.on.public.web) {
    .set.mSet(mSetObj)
    return(msg)
  } else {
    print(msg)
    return(.set.mSet(mSetObj))
  }
}

mSet <- PlotOPLS.Permutation(mSet, "OPLS-DA_permutation", 200, width = width, height = height)
file.rename("OPLS-DA_permutation.png", paste0(outPrefix, "_OPLS-DA_permutation.png"))
file.rename("OPLS-DA_permutation.pdf", paste0(outPrefix, "_OPLS-DA_permutation.pdf"))

mSet$analSet$oplsda$summaryDF %>%
  as_tibble %>%
  mutate(Group = str_glue("{basename(outPrefix)}")) %>%
  select(c(ncol(.), 2:ncol(.) - 1)) %>%
  select(-RMSEE) %>%
  write.xlsx(str_glue("{outPrefix}_OPLS-DA_summary.xlsx"))

print("Fig. 4B volcano plot of differential serum metabolites")

library(ggplot2)
library(tidyverse)

options(scipen = 100, stringsAsFactors = F)

infile <- "metabolite_volcano_input.txt"
rawdata <- read.table(infile, header = T, sep = "\t", quote = "", check.names = F)

plotdata <- rawdata %>% arrange(`P-value`)
FDR <- p.adjust(plotdata$`P-value`, method = "fdr", n = length(plotdata$`P-value`))
plotdata$FDR <- FDR

plotdata <- plotdata %>%
  mutate(Type = case_when(VIP > 1 & FDR < 0.05 & Log2FC > 0 ~ "up",
                          VIP > 1 & FDR < 0.05 & Log2FC < 0 ~ "down",
                          TRUE ~ "insig"))

num <- table(plotdata$Type)
insiglabel <- paste0("Nosig\n(n=", num["insig"], ")")
uplabel <- paste0("Up\n(n=", num["up"], ")")
downlabel <- paste0("Down\n(n=", num["down"], ")")
plotdata$Type[plotdata$Type == "insig"] <- insiglabel
plotdata$Type[plotdata$Type == "up"] <- uplabel
plotdata$Type[plotdata$Type == "down"] <- downlabel

p <- ggplot(plotdata, aes(VIP, -log10(`FDR`))) + geom_point(aes(color = Type)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "VIP scores", y = "-Log10(FDR)", color = " ") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(breaks = c(uplabel, insiglabel, downlabel), values = c("#a47e96", "grey", "#4c708f")) +
  theme_bw() + theme(text = element_text(size = 16),
                     title = element_text(size = 18),
                     axis.text = element_text(color = "black"),
                     legend.key.height = unit(1.2, "cm"),
                     panel.border = element_rect(size = 1.2))
cairo_pdf("volcano.pdf", width = 6.5, height = 5)
print(p)
dev.off()
png("volcano.png", width = 6.5, height = 5, units = "in", type = "cairo", res = 300)
print(p)
dev.off()

print("Fig. 4C KEGG pathway enrichment of differential metabolites at level 2")

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(openxlsx)

infile <- "KEGG_AddLevel_level2.xlsx"
data <- read_excel(infile, sheet = 1)

process_pathway_cids <- function(pathway_data) {
  all_cids <- character(0)

  for (i in 1:nrow(pathway_data)) {
    cid_string <- pathway_data$CIDList[i]

    if (is.na(cid_string) || cid_string == "") {
      next
    }

    cid_vector <- str_split(cid_string, pattern = "\\+")[[1]]

    cid_vector <- str_trim(cid_vector)

    all_cids <- c(all_cids, cid_vector[cid_vector != ""])
  }

  unique_cids <- unique(all_cids)

  return(sort(unique_cids))
}

result <- data %>%
  group_by(Kegg_level_2) %>%
  summarize(
    CID_count = {
      group_data <- cur_data()
      cids <- process_pathway_cids(group_data)
      length(cids)
    },
    CID_list = {
      group_data <- cur_data()
      cids <- process_pathway_cids(group_data)
      paste(cids, collapse = "+")
    }
  ) %>%
  ungroup()

result <- result %>%
  arrange(desc(CID_count))

print(head(result, 10))
print(paste0("Total pathway categories: ", nrow(result)))

max_cid_count <- max(result$CID_count)

detailed_df <- result

for (i in 1:max_cid_count) {
  col_name <- paste0("CID_", i)
  detailed_df[[col_name]] <- NA
}

for (i in 1:nrow(detailed_df)) {
  cid_list <- str_split(detailed_df$CID_list[i], pattern = "\\+")[[1]]

  for (j in 1:length(cid_list)) {
    col_name <- paste0("CID_", j)
    detailed_df[i, col_name] <- cid_list[j]
  }
}

cid_columns <- paste0("CID_", 1:max_cid_count)
detailed_df <- detailed_df[, c("Kegg_level_2", "CID_count", "CID_list", cid_columns)]

write.xlsx(result, "pathway_CID_summary.xlsx")

write.xlsx(detailed_df, "pathway_CID_summary_detailed.xlsx")

top_n <- 10

top10 <- result[1:top_n, ]

top10 <- top10 %>%
  arrange(CID_count) %>%
  mutate(Kegg_level_2 = factor(Kegg_level_2, levels = Kegg_level_2))

p <- ggplot(top10, aes(x = Kegg_level_2, y = CID_count)) +
  geom_bar(stat = "identity", fill = "gray", width = 0.7) +
  geom_text(aes(label = CID_count),
            hjust = -0.2,
            size = 3.5,
            color = "black") +
  coord_flip() +
  labs(
    title = "Number of DEMs",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(20, 40, 20, 20)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave("kegg_level2_bar.pdf", p, width = 6, height = 5)
ggsave("kegg_level2_bar.png", p, width = 6, height = 5, dpi = 300)

print("Fig. 4D co-inertia analysis of serum metabolites and clinical phenotypes")

library(ade4)
library(readxl)
library(ggplot2)
library(ggrepel)
library(latex2exp)
library(dplyr)
library(openxlsx)

infile <- "phenotype_metabolite.xlsx"
data <- read_excel(infile, col_types = "text", sheet = "metabolite")

group <- data$group

clinical_vars <- names(data)[3:151]

metabolite_start <- 152
metabolite_vars <- names(data)[metabolite_start:ncol(data)]

convert_to_numeric <- function(df, vars) {
  df_converted <- df
  for (var in vars) {
    df_converted[[var]] <- suppressWarnings(as.numeric(df[[var]]))
  }
  return(df_converted)
}

data_numeric <- convert_to_numeric(data, c(clinical_vars, metabolite_vars))

clinical_data <- data_numeric[, clinical_vars]
metabolite_data <- data_numeric[, metabolite_vars]

clinical_data <- clinical_data[, apply(clinical_data, 2, var, na.rm = TRUE) > 0]
metabolite_data <- metabolite_data[, apply(metabolite_data, 2, var, na.rm = TRUE) > 0]

complete_idx <- complete.cases(clinical_data, metabolite_data)
clinical_data <- clinical_data[complete_idx, ]
metabolite_data <- metabolite_data[complete_idx, ]
group_clean <- group[complete_idx]

print(paste0("Samples after filtering: ", nrow(clinical_data)))

clinical_scaled <- scale(clinical_data, center = TRUE, scale = TRUE)
metabolite_scaled <- scale(metabolite_data, center = TRUE, scale = TRUE)

pca_clinical <- dudi.pca(clinical_scaled, scannf = FALSE, nf = min(5, ncol(clinical_scaled)))
pca_metabolite <- dudi.pca(metabolite_scaled, scannf = FALSE, nf = min(5, ncol(metabolite_scaled)))

coin_all <- coinertia(pca_clinical, pca_metabolite, scannf = FALSE, nf = 2)
rand_all <- randtest(coin_all, nrepet = 99)
p_all <- rand_all$pvalue

print(paste0("Co-inertia P value: ", p_all))

mx <- coin_all$mX
my <- coin_all$mY

plot_data <- data.frame(
  mx1 = mx[, 1],
  mx2 = mx[, 2],
  my1 = my[, 1],
  my2 = my[, 2],
  group = group_clean
)

plot_data$group_label <- ifelse(plot_data$group == "1", "MD", "HC")
plot_data$group_label <- factor(plot_data$group_label, levels = c("MD", "HC"))

eig <- coin_all$eig
var_exp <- round(eig / sum(eig) * 100, 2)

p <- ggplot(plot_data) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray80", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray80", size = 0.5) +

  geom_segment(aes(x = 0, y = 0, xend = mx1, yend = mx2, color = group_label),
               alpha = 0.4, size = 0.3) +

  geom_segment(aes(x = 0, y = 0, xend = my1, yend = my2, color = group_label),
               alpha = 0.4, size = 0.3) +

  geom_point(aes(x = mx1, y = mx2, color = group_label, shape = "Clinical"),
             size = 4, alpha = 0.8) +

  geom_point(aes(x = my1, y = my2, color = group_label, shape = "Metabolite"),
             size = 4, alpha = 0.8) +

  stat_ellipse(aes(x = mx1, y = mx2, color = group_label),
               type = "norm", level = 0.95,
               linetype = "solid", size = 1, alpha = 0.3) +

  stat_ellipse(aes(x = my1, y = my2, color = group_label),
               type = "norm", level = 0.95,
               linetype = "solid", size = 1, alpha = 0.3) +

  scale_color_manual(
    values = c("MD" = "#E41A1C", "HC" = "#2C7BB6"),
    name = "Group"
  ) +

  scale_shape_manual(
    values = c("Clinical" = 16, "Metabolite" = 17),
    labels = c("Clinical" = "Clinical", "Metabolite" = "Metabolite"),
    name = "Variables"
  ) +

  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray95", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray60", size = 0.8),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "gray80", size = 0.5),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(size = 10, hjust = 1, color = "gray50")
  ) +

  xlab(paste0("Axis1 (", var_exp[1], "%)")) +
  ylab(paste0("Axis2 (", var_exp[2], "%)")) +

  labs(
    title = "Coia:MDD vs HC",
    subtitle = TeX(paste0("P = ", round(p_all, 4),
                          ", RV = ", round(coin_all$RV[1], 3),
                          ", n = ", nrow(plot_data))),
    caption = ""
  )

print(p)

ggsave("coia_serum_phenotype.pdf", p, width = 12, height = 10, dpi = 300)
ggsave("coia_serum_phenotype.png", p, width = 12, height = 10, dpi = 300)

write.xlsx(list(
  Scores = coin_all$li,
  Clinical_coords = coin_all$mX,
  Metabolite_coords = coin_all$mY,
  Eigenvalues = data.frame(Eigenvalue = coin_all$eig),
  Summary = data.frame(
    P_value = p_all,
    RV = coin_all$RV[1],
    N = nrow(plot_data),
    MD_n = sum(plot_data$group_label == "MD"),
    HC_n = sum(plot_data$group_label == "HC")
  )
), "coia_all_groups.xlsx")

print(paste0("Total samples: ", nrow(plot_data)))
print(paste0("MD samples: ", sum(plot_data$group_label == "MD")))
print(paste0("HC samples: ", sum(plot_data$group_label == "HC")))
print(paste0("P value: ", p_all))

print("Fig. 4E random forest and Spearman correlation of amino acids and derivatives")

library(psych)
library(reshape2)
library(ggplot2)
library(randomForest)
library(patchwork)
library(readxl)

infile <- "phenotype_amino_acid_correlation.xlsx"
data <- read_excel(infile, sheet = 1)

phenotype_cols <- c("HAMD", "PHQ", "HAMA", "GAD", "Spotter", "Symbol check",
                    "Codebreaker", "Trials", "PDQ5D", "TOC", "TCC",
                    "NEUT", "LYMPH", "MONO", "IL2", "IL4", "IL6", "IL10", "IL17", "IFNy", "TNFa")

phenotype_category <- data.frame(
  Phenotype = phenotype_cols,
  Category = c(
    rep("Emotional", 4),
    rep("Cognitive", 7),
    rep("Inflammatory", 10)
  ),
  stringsAsFactors = FALSE
)

amino_cols <- c("L-Glutamic acid", "Pyroglutamic acid", "N-Palmitoylglycine",
                "L-Threonine", "Taurine", "DL-Arginine", "N-Alpha-acetyllysine",
                "(±)-Tryptophan", "N-Phenylacetylphenylalanine", "NG-Amino-L-arginine",
                "5-hydroxy-L-tryptophan", "D-phenylalanine", "L-Histidinol",
                "N-Methylisoleucine", "L-Homophenylalanine", "N-Acetyl-L-tyrosine",
                "DL-Cysteine", "D-Pyroglutamic acid", "Amastatin",
                "O-benzyl-L-serine", "L-Formylkynurenine", "4-Hydroxy-l-tryptophan")

covariate_cols <- c("age", "BMI", "sex_1male")

phenotype_cols_exist <- phenotype_cols[phenotype_cols %in% colnames(data)]
amino_cols_exist <- amino_cols[amino_cols %in% colnames(data)]
covariate_cols_exist <- covariate_cols[covariate_cols %in% colnames(data)]

phenotype_category <- phenotype_category[phenotype_category$Phenotype %in% phenotype_cols_exist, ]

print(paste0("Phenotype columns: ", length(phenotype_cols_exist)))
print(paste0("Amino acid columns: ", length(amino_cols_exist)))

all_cols <- c(amino_cols_exist, phenotype_cols_exist, covariate_cols_exist)
new_names <- make.names(all_cols, unique = TRUE)
name_mapping <- setNames(new_names, all_cols)

data_renamed <- data[, all_cols, drop = FALSE]
colnames(data_renamed) <- new_names

amino_cols_clean <- name_mapping[amino_cols_exist]
phenotype_cols_clean <- name_mapping[phenotype_cols_exist]
covariate_cols_clean <- name_mapping[covariate_cols_exist]

phenotype_category$Phenotype_clean <- name_mapping[phenotype_category$Phenotype]

data_clean <- data_renamed
for (col in colnames(data_clean)) {
  data_clean[[col]] <- as.numeric(as.character(data_clean[[col]]))
}

data_clean <- na.omit(data_clean)
print(paste0("Rows after filtering: ", nrow(data_clean)))

cor_data <- data_clean[, c(amino_cols_clean, phenotype_cols_clean), drop = FALSE]

spearman_result <- corr.test(cor_data[, amino_cols_clean, drop = FALSE],
                             cor_data[, phenotype_cols_clean, drop = FALSE],
                             method = "spearman",
                             adjust = "none")

r <- data.frame(spearman_result$r)
r$Variable <- rownames(r)
r_melt <- melt(r, id = "Variable", variable.name = "Phenotype", value.name = "Correlation")

print(paste0("Correlation rows: ", nrow(r_melt)))

importance_list <- list()
r2_values <- c()
successful_pheno <- c()

for (i in 1:length(phenotype_cols_clean)) {
  pheno <- phenotype_cols_clean[i]
  pheno_orig <- names(name_mapping)[match(pheno, name_mapping)]
  print(paste0("Analysing: ", pheno_orig, " (", i, "/", length(phenotype_cols_clean), ")"))

  df <- data_clean[, c(amino_cols_clean, covariate_cols_clean, pheno), drop = FALSE]
  colnames(df)[ncol(df)] <- "y"

  if (!is.numeric(df$y)) {
    print("Warning: y is not numeric")
    next
  }

  if (sd(df$y, na.rm = TRUE) == 0) {
    print("Warning: y has no variation")
    next
  }

  tryCatch({
    set.seed(123)
    rf <- randomForest(y ~ ., data = df, importance = TRUE, ntree = 500)

    r2 <- rf$rsq[length(rf$rsq)]
    r2_values <- c(r2_values, r2)
    successful_pheno <- c(successful_pheno, pheno)

    print(paste0("R2 = ", round(r2, 4)))

    imp <- data.frame(importance(rf, scale = TRUE))
    imp$Variable <- rownames(imp)
    imp <- imp[imp$Variable %in% amino_cols_clean, ]

    if (nrow(imp) > 0) {
      imp_col <- grep("IncMSE", colnames(imp), value = TRUE)[1]
      if (is.na(imp_col)) imp_col <- colnames(imp)[1]

      imp$Phenotype <- pheno
      importance_list[[pheno]] <- imp[, c("Variable", "Phenotype", imp_col)]
      colnames(importance_list[[pheno]])[3] <- "Importance"
    }

  }, error = function(e) {
    print(paste0("Error: ", e$message))
  })
}

if (length(successful_pheno) == 0) {
  stop("No phenotype was analysed successfully, please check the input data")
}

print(paste0("Phenotypes analysed successfully: ", length(successful_pheno)))

importance_all <- do.call(rbind, importance_list)

r2_df <- data.frame(
  Phenotype = phenotype_cols_clean,
  R2 = NA
)

for (i in 1:length(successful_pheno)) {
  r2_df$R2[r2_df$Phenotype == successful_pheno[i]] <- r2_values[i]
}

r2_df <- r2_df[!is.na(r2_df$R2), ]

plot_data <- merge(r_melt, importance_all, by = c("Variable", "Phenotype"), all = TRUE)
plot_data$Importance[plot_data$Importance < 0] <- NA

plot_data$Variable_Display <- names(name_mapping)[match(plot_data$Variable, name_mapping)]
plot_data$Phenotype_Display <- names(name_mapping)[match(plot_data$Phenotype, name_mapping)]

pheno_order <- phenotype_cols_clean[phenotype_cols_clean %in% r2_df$Phenotype]
plot_data$Phenotype <- factor(plot_data$Phenotype, levels = pheno_order)
r2_df$Phenotype <- factor(r2_df$Phenotype, levels = pheno_order)

r2_df <- merge(r2_df, phenotype_category[, c("Phenotype_clean", "Category")],
               by.x = "Phenotype", by.y = "Phenotype_clean", all.x = TRUE)

var_order <- aggregate(Importance ~ Variable, data = plot_data, mean, na.rm = TRUE)
var_order <- var_order[order(var_order$Importance, decreasing = TRUE), "Variable"]
if (length(var_order) > 0) {
  plot_data$Variable <- factor(plot_data$Variable, levels = rev(var_order))
}

category_data <- data.frame(
  Phenotype = pheno_order,
  Category = r2_df$Category[match(pheno_order, r2_df$Phenotype)],
  y_position = max(as.numeric(plot_data$Variable)) + 0.5
)

category_colors <- c(
  "Emotional" = "#E41A1C",
  "Cognitive" = "#377EB8",
  "Inflammatory" = "#4DAF4A"
)

p_bar <- ggplot(r2_df, aes(x = Phenotype, y = R2 * 100)) +
  geom_col(fill = "steelblue", width = 0.7) +
  theme_bw() +
  theme_classic() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) + coord_cartesian(ylim = c(-10, 80)) +
  ylab("Explained variation(%)")

p_heatmap <- ggplot() +
  geom_tile(data = plot_data,
            aes(x = Phenotype, y = Variable, fill = Correlation)) +
  scale_fill_gradientn(
    colors = c("#2D6DB1", "white", "#DC1623"),
    limit = c(-1, 1),
    name = "Correlation"
  ) +
  geom_point(data = plot_data[!is.na(plot_data$Importance), ],
             aes(x = Phenotype, y = Variable, size = Importance),
             shape = 1) +
  scale_size_continuous(range = c(0, 6), name = "Importance (%)") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black"),
    legend.key = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 9),
    axis.ticks = element_line(color = "black"),
    axis.title = element_blank()
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

x_positions <- 1:length(pheno_order)
category_bar_data <- data.frame(
  x = x_positions,
  Phenotype = pheno_order,
  Category = r2_df$Category[match(pheno_order, r2_df$Phenotype)]
)

p_category <- ggplot(category_bar_data, aes(x = x, y = 1, fill = Category)) +
  geom_tile(height = 0.2) +
  scale_fill_manual(values = category_colors, name = "Category") +
  scale_x_continuous(
    breaks = x_positions,
    labels = names(name_mapping)[match(pheno_order, name_mapping)],
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               color = "black", size = 10),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  ) +
  labs(x = NULL, y = NULL)

combined_plot <- (p_bar / p_heatmap / p_category) +
  plot_layout(ncol = 1, heights = c(1, 5, 0.15)) +
  plot_annotation(
    title = "Random Forest Analysis with Spearman Correlation",
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
    )
  )

ggsave("RandomForest_Spearman_with_Category.pdf", combined_plot,
       width = 9, height = 9, dpi = 300)
ggsave("RandomForest_Spearman_with_Category.png", combined_plot,
       width = 9, height = 9, dpi = 300)

print(combined_plot)

write.csv(plot_data, "plot_data_final.csv", row.names = FALSE)
write.csv(r2_df, "r2_results_final.csv", row.names = FALSE)
