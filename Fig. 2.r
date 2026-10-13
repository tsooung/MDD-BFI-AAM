cl_list <- c("HC", "MDD")
for (cl in cl_list) {
  pl_list <- c("6")
  for (pl in pl_list) {
    cor_df = read.table(paste0(cl, "/", cl, "_bian", pl, ".csv"), header = T, sep = ",", comment.char = "", check.names = F)

    edgedata = read.table(paste0(cl, "/", cl, "_dian", pl, ".csv"), header = T, sep = ",", comment.char = "", check.names = F)

    library(igraph)
    library(WGCNA)

    igraph <- graph_from_data_frame(cor_df, directed = F, vertices = edgedata)

    nodes_num <- length(V(igraph))
    nodes_num

    edges_num <- length(E(igraph))
    edges_num

    average_degree <- mean(degree(igraph))
    average_degree <- 2 * edges_num / nodes_num
    average_degree

    nodes_connectivity <- vertex.connectivity(igraph)
    nodes_connectivity

    edges_connectivity <- edge.connectivity(igraph)
    edges_connectivity

    average_path_length <- average.path.length(igraph, directed = FALSE)
    average_path_length

    graph_diameter <- diameter(igraph, directed = FALSE)
    graph_diameter

    graph_density <- graph.density(igraph)
    graph_density

    clustering_coefficient <- transitivity(igraph)
    clustering_coefficient

    betweenness_centralization <- centralization.betweenness(igraph)$centralization
    betweenness_centralization

    degree_centralization <- centralization.degree(igraph)$centralization
    degree_centralization

    fc <- cluster_fast_greedy(igraph)
    modularity <- modularity(igraph, membership(fc))

    igraph_character <- data.frame(
      nodes_num,
      edges_num,
      average_degree,
      nodes_connectivity,
      edges_connectivity,
      average_path_length,
      graph_diameter,
      graph_density,
      clustering_coefficient,
      betweenness_centralization,
      degree_centralization,
      modularity
    )
    igraph_character

    newdir = paste0(cl, "/")

    write.table(igraph_character, paste0(newdir, "igraph_character.", cl, "_", pl, ".txt"), sep = '\t', row.names = FALSE, quote = FALSE)

  }
}

library(tidyverse)
library(circlize)

cl_list <- c("HC", "MDD")

pair_list <- list(
  c("Bacteria", "Eukaryota")
)

for (cl in cl_list) {
  print(cl)

  corr_df <- read.delim(paste0("HC MDD/median_correlation_", cl, ".tsv"), row.names = 1)
  pvalues_df <- read.delim(paste0("HC MDD/pvalues_", cl, ".tsv"), row.names = 1)

  corr_df <- rownames_to_column(corr_df, var = "Variable1")
  long_data_cor <- pivot_longer(corr_df, -Variable1, names_to = "Variable2", values_to = "Value") %>%
    filter(Variable1 < Variable2)

  pvalues_df <- rownames_to_column(pvalues_df, var = "Variable1")
  long_data_p <- pivot_longer(pvalues_df, -Variable1, names_to = "Variable2", values_to = "Value") %>%
    filter(Variable1 < Variable2)

  long_merge <- left_join(long_data_cor, long_data_p, by = c("Variable1", "Variable2"))
  colnames(long_merge) <- c("from", "to", "cor", "p")
  long_merge$padj <- p.adjust(long_merge$p, method = "BH")

  filtered <- long_merge %>%
    filter(padj <= 0.05, abs(cor) > 0.6)

  taxnomy <- read.table(
    "taxnomy/taxnomy_zong.txt",
    sep = "|",
    na.strings = c("", "NA"),
    stringsAsFactors = FALSE,
    header = FALSE
  ) %>%
    mutate(
      feature = gsub("s__", "", V8),
      kingdom = gsub("d__", "", V1)
    )

  df_anno <- filtered %>%
    left_join(taxnomy %>% select(feature, kingdom), by = c("from" = "feature")) %>%
    rename(from_kingdom = kingdom) %>%
    left_join(taxnomy %>% select(feature, kingdom), by = c("to" = "feature")) %>%
    rename(to_kingdom = kingdom)

  for (pair in pair_list) {
    k1 <- pair[1]
    k2 <- pair[2]

    bian <- df_anno %>%
      filter(
        (from_kingdom == k1 & to_kingdom == k2) |
          (from_kingdom == k2 & to_kingdom == k1)
      )

    out_dir <- file.path("HC MDD/", cl)
    write.csv(bian, file.path(out_dir, paste0(cl, "_bian6_", "Bacteria_Eukaryota", ".csv")),
              row.names = FALSE, quote = FALSE)
  }
}
