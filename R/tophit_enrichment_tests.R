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

Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)

base_folder <- here("results/full_run_ZengAWS.1698256033/") #2009 hits

merged_obs <- read_csv(file.path(base_folder, "merged.obs.csv.zip"))
universe_of_Zeng_ids <- merged_obs %>% filter(study_id == "Zeng") %>% pull(cell.type) %>% unique()
length(universe_of_Zeng_ids)
universe_of_Macosko_ids <- merged_obs %>% filter(study_id == "Macosko") %>% pull(cell.type) %>% unique()
length(universe_of_Macosko_ids)

data_base_folder <- here("data")
cluster_info <- read_tsv(here("data", "whole_mouse_brain", "zeng", "from_aws", "AIT21.0", "AIT21_annotation_freeze_081523.tsv" ))

cluster_info %>% pull(class_label) %>% unique() %>% length()
cluster_info %>% pull(subclass_label) %>% unique() %>% length()
cluster_info %>% group_by(class_label) %>% count() %>% as.data.frame()
cluster_info %>% group_by(subclass_label) %>% count() %>% as.data.frame() #too fine a resolution

high_level_size_threshold <- 1

cluster_info %<>% mutate(cl_with_study = paste0("Zeng|", cl ))

#remove LQ info
cluster_info %>% filter(!(cl %in% universe_of_Zeng_ids))
cluster_info %<>% filter(cl %in% universe_of_Zeng_ids)

#build GO sets for tmod -slow
tmodNames <- tibble()
highlevel_map <- list()

#create a tmod object with high level cell types for Zeng
for(high_level_name in cluster_info %>% pull(class_label) %>% unique()) {
  print(high_level_name)
  cluster_ids <- unique(cluster_info %>% filter(class_label == high_level_name) %>% pull(cl_with_study))  
  if(length(cluster_ids) < high_level_size_threshold) {  #filter for greater than five
    next
  }
  highlevel_map[high_level_name] <- list(cluster_ids)
  tmodNames <- bind_rows(tmodNames, data.frame(ID=high_level_name, Title = high_level_name))
}
Zeng_tmod_obj <- makeTmod(modules = tmodNames, modules2genes = highlevel_map)
Zeng_tmod_obj

Macosko_cluster_info_new <- read_csv(here("data", "whole_mouse_brain", "macosko", "from_google_drive", "CellType_MetaCluster.csv" ))
Macosko_cluster_info_new %<>% separate_rows(LeafLists, sep = "\\|")
Macosko_cluster_info_new %<>% mutate(Annotation = ifelse(is.na(Annotation), "", Annotation))
Macosko_cluster_info_new %<>% mutate(label = paste0(Meta_Cluster, " ", Annotation ))
Macosko_cluster_info_new %>% filter(Meta_Cluster == "MC_223")
Macosko_cluster_info_new %>% group_by(Meta_Cluster) %>% count() %>% as.data.frame()

Macosko_cluster_info_new %<>% select(label, annotation = LeafLists)

Macosko_cluster_info_new %>% group_by(label) %>% count() %>% as.data.frame()
colnames(Macosko_cluster_info_new)
Macosko_cluster_info_new %<>% mutate(annotation = paste0("Macosko|", annotation ))
high_level_size_threshold <- 10 #need at least many clusters
tmodNames <- tibble()
highlevel_map <- list()
#create a tmod object with high level cell types for Zeng
for(high_level_name in Macosko_cluster_info_new %>% pull(label) %>% unique()) {
  print(high_level_name)
  cluster_ids <- unique(Macosko_cluster_info_new %>% filter(label == high_level_name) %>% pull(annotation))  
  #cluster_ids <- paste0("Macosko|", cluster_ids) 
  if(length(cluster_ids) < high_level_size_threshold) {  #filter for greater than a threshold
    next
  }
  highlevel_map[high_level_name] <- list(cluster_ids)
  
  tmodNames <- bind_rows(tmodNames, data.frame(ID=high_level_name, Title = high_level_name))
}
Macosko_tmod_obj <- makeTmod(modules = tmodNames, modules2genes = highlevel_map)
Macosko_tmod_obj

write_out_tmod  <- function(results, base_folder, original_top_hits_filename, study, csv_suffix) {
  dir.create(paste0(base_folder, "top_hit_enrichments"), showWarnings = F)
  if ("E" %in% colnames(results)) {
    results %<>%
      mutate(
        expected = n * B / N,
        bias   = case_when(
          b > expected ~ "enrichment",
          b < expected ~ "depletion")) %>% select(-expected) 
    results %<>% mutate(
      # Two-sided Fisher exact p-value (sum of probs <= observed table), exact test
      two_sided.P.Value = mapply(function(b, B, n, N) {
        # 2x2 table:
        #                In set   Not in set
        # In sample         b       n - b
        # Not in sample   B - b   N - n - B + b
        mat <- matrix(c(b, n - b, B - b, N - n - B + b), nrow = 2, byrow = TRUE)
        fisher.test(mat, alternative = "two.sided")$p.value
      }, b, B, n, N))
    results %<>% select(-ID, -E)
    results %<>% rename(Clusters = B)
    results %<>% rename(Overlap = b)
    results %<>% mutate(one_sided.P.Value = signif(P.Value, digits = 2), one_sided.adj.P.Val = signif(adj.P.Val, digits = 2),
                        two_sided.P.Value = signif(two_sided.P.Value, digits = 2), two_sided.adj.P.Val = signif(p.adjust(two_sided.P.Value,method = 'fdr'), digits = 2)
    ) %>% select(-P.Value, -adj.P.Val)
    results %<>% select(-one_sided.P.Value, -one_sided.adj.P.Val) #just go with two sided tests in the end - could be refactored to remove tmod
    #results %<>% relocate(two_sided.P.Value, two_sided.adj.P.Val, .after = last_col())
    
  } else if ("AUC" %in% colnames(results)) {
    results %<>% select(-ID, -U)
    results %<>% mutate(P.Value = P.Value *2 ) # double because tmod runs one sided tests
    results %<>% mutate(adj.P.Val = p.adjust(P.Value, method = "fdr") )
    results %<>% mutate(P.Value = signif(P.Value, digits = 2), adj.P.Val = signif(adj.P.Val, digits = 2), AUC = signif(AUC, digits = 2))
    results %<>% rename(AUROC = AUC, Count=N1)
  }
  results %>% write_csv(paste0(base_folder, "top_hit_enrichments/", gsub(".csv", paste0(".", study, csv_suffix), original_top_hits_filename)))
}

all_recip_results <- tibble()
all_non_matching_results <- tibble()

#top_hit_files <- list.files(path = base_folder, pattern = "op_hit") #all possible
top_hit_files <- c("top_hits.0.95.csv") #just run one for simplicity

for (top_hit_filename in top_hit_files) {
  if (top_hit_filename == "top_hit_enrichments") {
    next
  }
  if (grepl("with_sizes", top_hit_filename)) {
    next
  }
  top_hits <- read_csv(paste0(base_folder, top_hit_filename))
  if ("...1" %in% colnames(top_hits))  {
    top_hits %<>% select(-`...1`)
  }
  
  top_hits %<>% mutate(Celltype_1_study = getStudyId(`Study_ID|Celltype_1`))
  top_hits %<>% mutate(Celltype_2_study = getStudyId(`Study_ID|Celltype_2`))
  all_cell_types_with_hits <- union(top_hits$`Study_ID|Celltype_1`, top_hits$`Study_ID|Celltype_2`)
  
  #with recip hits
  all_recip_hits <- top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% select(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`) %>% pivot_longer(everything()) %>% pull(value)
  length(all_recip_hits)
  result <- NULL
  result <- tmodHGtest(fg = all_recip_hits, bg = Zeng_tmod_obj$gv, mset = Zeng_tmod_obj, qval = 1.01, filter = TRUE) %>% tibble()
  write_out_tmod(result, base_folder, top_hit_filename, "Zeng", ".reciprocal_hits_enrichment.csv")
  all_recip_results %<>% bind_rows(result %>% mutate(study = "Zeng", filename = top_hit_filename))
  
  result <- NULL
  result <- tmodHGtest(fg = all_recip_hits, bg = Macosko_tmod_obj$gv, mset = Macosko_tmod_obj, qval = 1.01, filter = TRUE) %>% tibble()
  write_out_tmod(result, base_folder, top_hit_filename, "Macosko", ".reciprocal_hits_enrichment.csv")
  result %<>% mutate(Title = gsub("^MC_[0-9]+ ", "", Title))
  write_out_tmod(result, base_folder, top_hit_filename, "Macosko_MC_removed", ".reciprocal_hits_enrichment.csv")
  all_recip_results %<>% bind_rows(result %>% mutate(study = "Macosko", filename = top_hit_filename))
  
  
  
  #Macosko types without tophits
  Macosko_without <- setdiff(Macosko_tmod_obj$gv, all_cell_types_with_hits)
  
  
  #Zeng types without tophits
  Zeng_without <- setdiff(Zeng_tmod_obj$gv, all_cell_types_with_hits)
  length(Zeng_without)
  #annotation
  
}

#stats on those without hits
merged_obs %>% group_by(study_id, cell.type) %>% count() %>% group_by(study_id) %>% summarize(median = median(n))

merged_obs %>% group_by(study_id, cell.type) %>% count() %>% mutate(combined = paste0(study_id, "|", cell.type)) %>% filter(combined %in% Macosko_without)
merged_obs %>% group_by(study_id, cell.type) %>% count() %>% mutate(combined = paste0(study_id, "|", cell.type)) %>% filter(combined %in% Macosko_without)%>% group_by(study_id) %>% summarize(median = median(n))


merged_obs %>% group_by(study_id, cell.type) %>% count() %>% mutate(combined = paste0(study_id, "|", cell.type)) %>% filter(combined %in% Zeng_without) %>% as.data.frame() %>% arrange(-n)
merged_obs %>% group_by(study_id, cell.type) %>% count() %>% mutate(combined = paste0(study_id, "|", cell.type)) %>% filter(combined %in% Zeng_without) %>% as.data.frame() %>% arrange(-n) %>% summarize(median = median(n))

