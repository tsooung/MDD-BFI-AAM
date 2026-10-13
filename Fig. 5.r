library(ade4)
library(readxl)

data <- read_excel("paper5.xlsx", col_types = "text")

convert_to_numeric <- function(df) {
  df_converted <- df
  for (i in 1:ncol(df)) {
    if (is.character(df[[i]])) {
      numeric_test <- suppressWarnings(as.numeric(df[[i]]))
      if (!all(is.na(numeric_test))) {
        df_converted[[i]] <- numeric_test
      }
    }
  }
  return(df_converted)
}

data_numeric <- convert_to_numeric(data)

n_env <- 14
X <- data_numeric[, 4:n_env]
Y <- data_numeric[, (n_env + 1):50]

complete_cases <- complete.cases(X, Y)
X <- X[complete_cases, ]
Y <- Y[complete_cases, ]

pca_X <- dudi.pca(X, scannf = FALSE, nf = 3)
pca_Y <- dudi.pca(Y, scannf = FALSE, nf = 3)

coin_result <- coinertia(pca_X, pca_Y, scannf = FALSE, nf = 2)

summary(coin_result)

eigenvalues <- coin_result$eig
axes <- coin_result$co
scores <- coin_result$li

plot(coin_result)
par(mfrow = c(1, 2))
plot(coin_result, xax = 1, yax = 2)

required_packages <- c(
  "readxl", "dplyr", "ggplot2", "ggalluvial", "ggnewscale"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  install.packages(missing_packages, dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  library(ggnewscale)
})

input_file <- "mediation_results_final.xlsx"
result_sheet <- "plot"

output_dir <- "final_figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

microbe_map <- data.frame(
  Bug = c(
    "Zoopagomycota", "Oomycota", "Microsporidia", "Blastocladiomycota", "Chytridiomycota",
    "Fibrobacterota", "Planctomycetota", "Gemmatimonadota", "Campylobacterota", "Cyanobacteriota"
  ),
  Microbe_type = c(rep("fungi", 5), rep("bacteria", 5)),
  stringsAsFactors = FALSE
)

phenotype_group_map <- data.frame(
  Phenotype = c(
    "HAMD_Total", "PHQ9_Total", "HAMA_Total", "GAD7_Total",
    "Spotter", "Symbol_Check", "Codebreaker", "Trials", "PDQ_5_D", "THINC_Obj", "THINC_Comp",
    "MONO_Pct", "IL_4", "IL_10", "IL_17", "TNF_alpha"
  ),
  Phenotype_group = c(
    rep("Emotion", 4),
    rep("Cognition", 7),
    rep("Immune", 5)
  ),
  stringsAsFactors = FALSE
)

phenotype_order <- phenotype_group_map$Phenotype

results_raw <- readxl::read_xlsx(input_file, sheet = result_sheet)
required_columns <- c("Bug", "Metabolite", "Phenotype",
                      "ACME_std")
missing_columns <- setdiff(required_columns, names(results_raw))
if (length(missing_columns) > 0L) {
  stop("The result table is missing the following columns: ", paste(missing_columns, collapse = ", "))
}

confirmed_results <- results_raw |>
  dplyr::select(Bug, Metabolite, Phenotype, dplyr::everything()) |>
  dplyr::left_join(microbe_map, by = "Bug") |>
  dplyr::left_join(phenotype_group_map, by = "Phenotype") |>
  dplyr::mutate(
    Effect_direction = dplyr::case_when(
      ACME_std > 0 ~ "positive",
      ACME_std < 0 ~ "negative",
      TRUE ~ "zero"
    )
  )

if (nrow(confirmed_results) == 0L) stop("No mediation path was read from the result table.")

unmapped_bugs <- confirmed_results |>
  dplyr::filter(is.na(Microbe_type)) |>
  dplyr::distinct(Bug) |> dplyr::pull(Bug)
if (length(unmapped_bugs) > 0L) {
  stop("The following phyla have no bacteria/fungi label, please extend microbe_map: ",
       paste(unmapped_bugs, collapse = ", "))
}

unmapped_pheno <- confirmed_results |>
  dplyr::filter(is.na(Phenotype_group)) |>
  dplyr::distinct(Phenotype) |> dplyr::pull(Phenotype)
if (length(unmapped_pheno) > 0L) {
  stop("The following phenotypes have no group label, please extend phenotype_group_map: ",
       paste(unmapped_pheno, collapse = ", "))
}

sankey_data <- confirmed_results |>
  dplyr::count(Microbe_type, Effect_direction, Bug, Metabolite, Phenotype, name = "n_paths")

bug_order <- c(
  "Fibrobacterota", "Planctomycetota", "Gemmatimonadota", "Campylobacterota", "Cyanobacteriota",
  "Zoopagomycota", "Oomycota", "Microsporidia", "Blastocladiomycota", "Chytridiomycota"
)

sankey_data <- sankey_data |>
  dplyr::mutate(
    Bug = factor(Bug, levels = bug_order),
    Microbe_type = factor(Microbe_type, levels = c("bacteria", "fungi")),
    Effect_direction = factor(Effect_direction, levels = c("positive", "negative", "zero")),
    Phenotype = factor(Phenotype, levels = phenotype_order)
  )

bacteria_color <- "#75C2AE"
fungi_color    <- "#FFE090"
positive_color <- "#EE7371"
negative_color <- "#7190CB"
zero_color     <- "#9E9E9E"

emotion_color  <- "#D4B8E8"
cognition_color <- "#A8D8D0"
immune_color   <- "#C9D3E0"
metabolite_color <- "#E9EDF2"

node_color_values <- c(
  bacteria = bacteria_color,
  fungi = fungi_color,
  Emotion  = emotion_color,
  Cognition = cognition_color,
  Immune   = immune_color
)

for (bug_name in unique(as.character(sankey_data$Bug))) {
  if (is.na(bug_name)) next
  bug_type <- as.character(microbe_map$Microbe_type[match(bug_name, microbe_map$Bug)])
  node_color_values[bug_name] <- if (identical(bug_type, "bacteria")) bacteria_color else fungi_color
}

for (pheno_name in unique(as.character(sankey_data$Phenotype))) {
  if (is.na(pheno_name)) next
  grp <- as.character(phenotype_group_map$Phenotype_group[match(pheno_name, phenotype_group_map$Phenotype)])
  node_color_values[pheno_name] <- switch(
    grp,
    Emotion = emotion_color,
    Cognition = cognition_color,
    Immune = immune_color,
    metabolite_color
  )
}

for (met_name in unique(as.character(sankey_data$Metabolite))) {
  if (is.na(met_name)) next
  node_color_values[met_name] <- metabolite_color
}

bacteria_legend_key <- sankey_data |>
  dplyr::filter(Microbe_type == "bacteria") |>
  dplyr::pull(Bug) |> as.character() |> unique() |> stats::na.omit() |> head(1)
fungi_legend_key <- sankey_data |>
  dplyr::filter(Microbe_type == "fungi") |>
  dplyr::pull(Bug) |> as.character() |> unique() |> stats::na.omit() |> head(1)

emotion_legend_key <- sankey_data |>
  dplyr::left_join(phenotype_group_map, by = "Phenotype") |>
  dplyr::filter(Phenotype_group == "Emotion") |>
  dplyr::pull(Phenotype) |> as.character() |> unique() |> stats::na.omit() |> head(1)
cognition_legend_key <- sankey_data |>
  dplyr::left_join(phenotype_group_map, by = "Phenotype") |>
  dplyr::filter(Phenotype_group == "Cognition") |>
  dplyr::pull(Phenotype) |> as.character() |> unique() |> stats::na.omit() |> head(1)
immune_legend_key <- sankey_data |>
  dplyr::left_join(phenotype_group_map, by = "Phenotype") |>
  dplyr::filter(Phenotype_group == "Immune") |>
  dplyr::pull(Phenotype) |> as.character() |> unique() |> stats::na.omit() |> head(1)

node_legend_breaks <- c(bacteria_legend_key, fungi_legend_key,
                        emotion_legend_key, cognition_legend_key, immune_legend_key)
node_legend_labels <- c("Bacteria", "Fungi", "Emotion", "Cognition", "Immune")

p <- ggplot2::ggplot(
  sankey_data,
  ggplot2::aes(axis1 = Bug, axis2 = Metabolite, axis3 = Phenotype, y = n_paths)
) +
  ggalluvial::geom_alluvium(
    ggplot2::aes(fill = Effect_direction),
    width = 0.22, alpha = 0.62, knot.pos = 0.35, color = NA
  ) +
  ggplot2::scale_fill_manual(
    values = c(positive = positive_color, negative = negative_color, zero = zero_color),
    breaks = c("positive", "negative"),
    name = "ACME", drop = FALSE,
    guide = ggplot2::guide_legend(override.aes = list(alpha = 0.85), order = 1)
  ) +
  ggnewscale::new_scale_fill() +
  ggalluvial::geom_stratum(
    ggplot2::aes(fill = ggplot2::after_stat(stratum)),
    width = 0.22, color = "white", linewidth = 0.45
  ) +
  ggplot2::geom_text(
    stat = "stratum",
    ggplot2::aes(label = ggplot2::after_stat(stratum)),
    size = 3.0, color = "#202124", lineheight = 0.95
  ) +
  ggplot2::scale_x_discrete(
    limits = c("Domain (X)", "Metabolite (M)", "Phenotype (Y)"),
    expand = c(0.10, 0.10)
  ) +
  ggplot2::scale_fill_manual(
    values = node_color_values,
    breaks = node_legend_breaks,
    labels = node_legend_labels,
    name = "domain", drop = FALSE,
    guide = ggplot2::guide_legend(order = 2, nrow = 1)
  ) +
  ggplot2::theme_void(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 17, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 6)),
    plot.subtitle = ggplot2::element_text(size = 11, color = "#4D5156", hjust = 0.5, margin = ggplot2::margin(b = 12)),
    axis.text.x = ggplot2::element_text(size = 12, face = "bold", color = "#202124", margin = ggplot2::margin(t = 8)),
    legend.position = "top",
    legend.justification = "center",
    legend.title = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(size = 9, color = "#6B7280", hjust = 0),
    plot.margin = ggplot2::margin(18, 35, 18, 35)
  )

png_file <- file.path(output_dir, "Sankey_98_Bug_Metabolite_Phenotype_ACME_direction.png")
pdf_file <- file.path(output_dir, "Sankey_98_Bug_Metabolite_Phenotype_ACME_direction.pdf")
csv_file <- file.path(output_dir, "Sankey_98_path_counts_with_direction.csv")

ggplot2::ggsave(filename = png_file, plot = p, width = 8, height = 11, units = "in", dpi = 600, bg = "white")
ggplot2::ggsave(filename = pdf_file, plot = p, width = 8, height = 11, units = "in", device = grDevices::cairo_pdf, bg = "white")
utils::write.csv(sankey_data, csv_file, row.names = FALSE, fileEncoding = "UTF-8")

print(p)

message("======================================================")
message("Three-column Sankey diagram finished")
message("Number of confirmed paths: ", sum(sankey_data$n_paths))
message("Bacteria-initiated paths: ", sum(sankey_data$n_paths[sankey_data$Microbe_type == "bacteria"]))
message("Fungi-initiated paths: ", sum(sankey_data$n_paths[sankey_data$Microbe_type == "fungi"]))
message("Positive ACME paths: ", sum(sankey_data$n_paths[sankey_data$Effect_direction == "positive"]))
message("Negative ACME paths: ", sum(sankey_data$n_paths[sankey_data$Effect_direction == "negative"]))
message("PNG: ", png_file)
message("PDF: ", pdf_file)
message("Path table: ", csv_file)
message("======================================================")
