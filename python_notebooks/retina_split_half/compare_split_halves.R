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

#preprint data
#base_folder <- "/vault/lfrench/temp_transfer/retina_cellxgene/split_half.1774982548/"
#published data
base_folder <- "/vault/lfrench/temp_transfer/retina_cellxgene/split_half.1775145180/"


top_hits <- read_csv(paste0(base_folder, "/top_hits.0.95.csv"))
recip_types <- union(top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_1`),
                     top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_2`))

#Count HVGs
#count cells
merged_obs <- read_csv(paste0(base_folder, "/merged.obs.csv.zip"))
merged_var <- read_csv(paste0(base_folder, "/merged.var.csv.zip"))
# nrow(merged_var)

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


top_hits %<>% mutate(is_same_type = getCellType(`Study_ID|Celltype_1`) == getCellType(`Study_ID|Celltype_2`)) 

## get the overlap of matched clusters that participate in the 2009
matched_particpaing_clusters <- top_hits %>% filter(Match_type == "Reciprocal_top_hit", is_same_type)
nrow(matched_particpaing_clusters)
top_hits %>% arrange(Mean_AUROC)

#get sizes of matched clusters
for_sizes <- top_hits
cluster_sizes <- merged_obs %>% group_by(paste0(study_id,"|", cell.type)) %>% count() %>% rename(cell.type = `paste0(study_id, "|", cell.type)`)
cluster_sizes %>% arrange(n)
cluster_sizes %<>% mutate(is_matched = cell.type %in% (matched_particpaing_clusters %>% pull(`Study_ID|Celltype_1`)) | cell.type %in% (matched_particpaing_clusters %>% pull(`Study_ID|Celltype_2`)))
cluster_sizes %>% group_by(is_matched) %>% summarize(median_size = median(n))
cluster_sizes %>% group_by(is_matched) %>% summarize(n=n())

median_size_unmatched <- cluster_sizes %>% filter(!is_matched) %>% pull(n) %>% median()
cluster_sizes %>% filter(!is_matched) %>% pull(n) %>% mean()
median_size_matched   <- cluster_sizes %>% filter(is_matched) %>% pull(n) %>% median()
cluster_sizes %>% filter(is_matched) %>% pull(n) %>% mean()
cluster_sizes %>% arrange(n)

