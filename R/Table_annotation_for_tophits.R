library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(here)
library(magrittr)
library(conflicted)
conflict_prefer("first", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")

#load top hits, filter, add cell counts, add names
base_folder <- here("results/full_run_ZengAWS.1698256033/") 

#uncomment/comment the below to switch tables
#top_hits <- read_csv(paste0(base_folder, "top_hits.0.95.csv"))
#out_filename <- "Supplement_table.top_recip_hits.csv"

top_hits <- read_csv(paste0(base_folder, "top_hits_asymmetric.0.99.best_vs_next.0.6.filtered.csv")) #612 high confidence
out_filename <- "Supplement_table.top_recip_hits.612.csv"

top_hits %<>% filter(Match_type == "Reciprocal_top_hit")
top_hits %>% nrow()

top_hits %<>% select(-`...1`)
top_hits %<>% rowwise() %>% mutate(first = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[2],
                                   second = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[1])
top_hits %<>% rename(cl_singlecell = first, name_singlenucleus = second)
top_hits %<>% select(-`Study_ID|Celltype_1`, -`Study_ID|Celltype_2`)
top_hits

merged_obs <- read_csv(paste0( base_folder, "merged.obs.csv.zip"))
merged_obs %<>% rename(cell_id = "...1")
cell_type_sizes <- merged_obs %>% group_by(cell.type, study_id) %>% count()
cell_type_sizes %<>% ungroup()
cell_type_sizes %<>% mutate(cell.type = paste0(study_id, "|", cell.type)) %>% select(-study_id)

#merge with tophits
top_hits %<>% inner_join(cell_type_sizes %>% select(cl_singlecell = cell.type, cell_count_singlecell = n))
top_hits %<>% inner_join(cell_type_sizes %>% select(name_singlenucleus = cell.type, cell_count_singlenucleus = n))

#add in cluster_id
Zeng_cluster_info <- read_tsv(here("data/whole_mouse_brain/zeng/from_aws/AIT21.0/AIT21_annotation_freeze_081523.tsv"))
Zeng_cluster_info %<>% mutate(cl = as.character(cl))
Zeng_cluster_info %<>% select(cl_singlecell = cl, cluster_id_singlecell = cluster_id, cluster_label_singlecell = cluster_id_label)
Zeng_cluster_info %<>% mutate(cluster_id_singlecell = as.character(cluster_id_singlecell))
Zeng_cluster_info %<>% mutate(cluster_id_singlecell = str_pad(cluster_id_singlecell, width = 4, side = "left", pad = "0"))
Zeng_cluster_info %<>% mutate(knowledgeID_singlecell = paste0("CS20230722_CLUS_", cluster_id_singlecell))

#remove prefixes
top_hits %<>% mutate(cl_singlecell = gsub(".*[|]", "", cl_singlecell))
top_hits %<>% mutate(name_singlenucleus = gsub(".*[|]", "", name_singlenucleus))

top_hits %<>% inner_join(Zeng_cluster_info)
top_hits[1,] %>% as.data.frame()

colnames(top_hits)
top_hits %<>% select(cl_singlecell, cluster_id_singlecell, knowledgeID_singlecell, cluster_label_singlecell, 
                     name_singlenucleus, cell_count_singlecell, cell_count_singlenucleus, Mean_AUROC)

top_hits %>% write_csv(paste0(base_folder, out_filename))
