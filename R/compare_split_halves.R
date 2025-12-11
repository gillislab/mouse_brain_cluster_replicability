library(ggExtra)
library(ggplot2)
library(tidyr)
library(MetaNeighbor)
library(conflicted)
library(magrittr)
library(reshape2)
library(dplyr)
library(here)
library(readr)
library(tmod)
conflict_prefer("count", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflicts_prefer(dplyr::rename)
conflicts_prefer(base::union)
conflicts_prefer(dplyr::first)
#count clusters - good check


base_folder_original <- here("results/full_run_ZengAWS.1698256033/") #2009 reciprocal result
top_hits <- read_csv(paste0(base_folder_original, "/top_hits.0.95.csv")) %>% filter(Match_type == "Reciprocal_top_hit")
#get the original top hits as strings
pair_strings_original <- top_hits %>% filter(Match_type == "Reciprocal_top_hit")  %>%
  rowwise() %>%
  mutate(
    pair = paste(
      sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`)),
      collapse = " || "
    )
  ) %>%
  ungroup() %>%
  pull(pair)


dirs <- list.dirs(path = here("results"), full.names = TRUE, recursive = FALSE)
dirs <- dirs[grepl("^split_half", basename(dirs))]
dirs
all_results <- tibble()

for(base_folder in dirs) {
  ##Count recip hits
  top_hits <- read_csv(paste0(base_folder, "/top_hits.0.95.csv"))
  recip_types <- union(top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_1`),
                       top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_2`))
  
  #Count HVGs
  #count cells
  merged_obs <- read_csv(paste0(base_folder, "/merged.obs.csv.gz"))
  merged_var <- read_csv(paste0(base_folder, "/merged.var.csv.gz"))
  # nrow(merged_var)
  # nrow(merged_obs)

  # # Get cell type sizes
  # cell_type_sizes <- merged_obs %>%
  #   group_by(cell.type = paste0(study_id, "|", cell.type)) %>%
  #   count()
  
  # merged_obs %>%
  #   group_by(study_id) %>%
  #   count()
  # 
  # merged_obs %>% select(study_id, cell.type) %>% distinct() %>%
  #   group_by(study_id) %>%
  #   count()
  # 
  wider_counts <- merged_obs %>%
    distinct(study_id, cell.type) %>%
    count(study_id, name = "n") %>%
    # optionally sort so "first" and "second" are alphabetical:
    # arrange(study_id) %>%
    summarise(
      study_id_1 = first(study_id),
      study_id_2 = nth(study_id, 2),
      n_clusters_1        = first(n),
      n_clusters_2        = nth(n, 2)
    )

  #get the top hits as strings
  pair_strings <- top_hits %>% filter(Match_type == "Reciprocal_top_hit")  %>%
    rowwise() %>%
    mutate(
      pair = paste(
        sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`)),
        collapse = " || "
      )
    ) %>%
    ungroup() %>%
    pull(pair)
  pair_strings <- gsub("_[BA]\\|", "|", pair_strings)
  
  single_result <- tibble(recip_hits = length(recip_types)/2, recip_hits_tophits = top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% nrow(), recip_pair_strings = length(pair_strings), cells = nrow(merged_obs), HVGs = nrow(merged_var),
                          matching_cluster_recip_hits = top_hits %>% mutate(is_same_type = getCellType(`Study_ID|Celltype_1`) == getCellType(`Study_ID|Celltype_2`)) %>% filter(Match_type == "Reciprocal_top_hit", is_same_type) %>% nrow(),
                          overlap_with_2009 = length(intersect(pair_strings_original, pair_strings)))
  single_result <- bind_cols(wider_counts, single_result)
  all_results %<>% bind_rows(single_result)
}
all_results
all_results %<>% mutate(overlap_with_2009_prop = overlap_with_2009/ length(pair_strings_original))
all_results %>% filter(n_clusters_1 != n_clusters_2) %>% pull(recip_hits) %>% mean()
all_results %>% filter(n_clusters_1 != n_clusters_2) %>% pull(recip_hits) %>% median()

all_results %>% filter(n_clusters_1 != n_clusters_2) %>% pull(overlap_with_2009_prop) %>% mean()



#plot

# ---- load packages ----
library(tidyverse)   # dplyr + tibble + readr + ggplot2, etc.
library(tidygraph)   # tidy interface to igraph objects
library(ggraph)      # grammar of graphics for graphs


edges <- all_results


# ---- build a directed tbl_graph ----
g <- as_tbl_graph(edges, directed = TRUE)

# ---- plot with ggraph ----
ggraph(g, layout = "graphopt") +
  # draw edges, mapping 'recip_hits' to both width and label
  geom_edge_link(aes(width = recip_hits,
                     label = paste0("                 ",recip_hits)),
                 angle_calc  = "along",
                 label_dodge = unit(2.5, "mm"),
                 #arrow       = arrow(length = unit(4, "mm")),
                 #end_cap     = circle(3, "mm"),
                 show.legend = FALSE) +
  scale_edge_width(range = c(0.3, 2)) +   # thinner to thicker edges
  # draw nodes
  geom_node_point(size = 5) +
  geom_node_text(aes(label = name), 
                 vjust = -1.2, size = 4) +
  theme_graph() +                         # clean theme for network plots
  labs(title = "recip_hits between single-cell studies")
