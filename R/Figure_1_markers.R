library(zoo)
library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(readxl)
library(here)
library(magrittr)
library(MetaNeighbor)
library(reshape2)
library(gplots)
library(conflicted)
conflict_prefer("first", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("union", "base")
conflict_prefer("setdiff", "base")

Sys.setenv("VROOM_CONNECTION_SIZE"= 21000000)
Macosko_results <- read_csv(here("results", "marker_results", "markers_on_Macosko.1712974391.csv.gz"))
Macosko_results %<>% rename(cluster_name = marker_name) 
colnames(Macosko_results) %>% head()

####
# load tophits
####
base_folder <- here("results", "full_run_ZengAWS.1698256033/") #updated Allen data - 2009 recips

top_hits <- read_csv(paste0(base_folder, "top_hits.0.95.csv"))

#check tophits length
c(top_hits %>% pull(`Study_ID|Celltype_1`), top_hits %>% pull(`Study_ID|Celltype_2`)) %>% unique() %>% length()

top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% nrow()

top_hits %<>% select(-`...1`)
top_hits %<>% rowwise() %>% mutate(first = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[2],
                                   second = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[1])
top_hits %<>% mutate(`Study_ID|Celltype_1` = first, `Study_ID|Celltype_2` = second)
top_hits %<>% select(-first, -second)
top_hits %>% group_by(Match_type) %>% count()
top_hits %>% nrow()

recip_types <- union(top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_1`) %>% as.character(),
                     top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_2`) %>% as.character())
#############
#############

###################
## Melt Macosko
###################
colnames(Macosko_results) %>% head()

Macosko_results_melted <- melt(Macosko_results) %>% tibble()
Macosko_results_melted %<>% rename(marker_set = variable)
Macosko_results_melted %<>% rename(auroc = value)
Macosko_results_melted

#Some code duplication below
best_clusters_per_marker_Mac <- Macosko_results_melted %>% group_by(marker_set) %>% arrange(-auroc) %>% summarize(best_cluster = first(cluster_name), best_AUROC = first(auroc),
                                                                                                                  second_cluster = nth(cluster_name, 2), second_AUROC = nth(auroc,2))
best_clusters_per_marker_Mac %<>% mutate(marker_source = if_else(grepl("^[0-9]+_markers$", marker_set), "Zeng", "Macosko"))
best_clusters_per_marker_Mac %<>% mutate(best_cluster_source = if_else(grepl("^[0-9]+_cluster$", best_cluster), "Zeng", "Macosko"))
best_clusters_per_marker_Mac

Zeng_markers_on_Macosko <- best_clusters_per_marker_Mac %>% filter(marker_source == "Zeng", best_cluster_source == "Macosko")
Zeng_markers_on_Macosko %>% pull(best_cluster) %>% unique() %>% length()
Zeng_markers_on_Macosko %>% pull(marker_set) %>% unique() %>% length()

#one to one mappings based on Zeng markers, can be used for spatial testing
Zeng_markers_on_Macosko_one_to_one <- Zeng_markers_on_Macosko %>% group_by(best_cluster) %>% slice_max(best_AUROC, with_ties = FALSE) 
Zeng_markers_on_Macosko_one_to_one
Zeng_markers_on_Macosko_one_to_one %<>% mutate(marker_set = gsub("_markers", "", marker_set))
Zeng_markers_on_Macosko_one_to_one %<>% mutate(best_cluster = gsub("_cluster", "", best_cluster))
Zeng_markers_on_Macosko_one_to_one %<>% select(marker_set, best_cluster, best_AUROC)

Zeng_markers_on_Macosko_one_to_one %>% write_csv(here("results", "marker_results", "Zeng_markers_on_Macosko_one_to_one_mapping.csv"))
marker_hits <- Zeng_markers_on_Macosko_one_to_one
marker_hits %<>% rename(`Study_ID|Celltype_1` = marker_set, `Study_ID|Celltype_2` = best_cluster)
marker_hits %<>% mutate(`Study_ID|Celltype_1` = paste0("Zeng|", as.character(`Study_ID|Celltype_1`)))
marker_hits %<>% mutate(`Study_ID|Celltype_2` = paste0("Macosko|", as.character(`Study_ID|Celltype_2`)))
#write out as a decoy tophits file for testing
marker_hits %>% mutate(Mean_AUROC = best_AUROC, Match_type = "Reciprocal_top_hit") %>% 
  write_csv(here("results/marker_results/Top_hit_format_for_Zeng_markers_on_Macosko_one_to_one_mapping.csv"))


nrow(Zeng_markers_on_Macosko_one_to_one)

####################################
#using full Zeng universe of genes:
####################################
Zeng_results <- read_csv(here("results", "marker_results", "markers_on_Zeng.1713323744.csv.gz"))

Zeng_results_melted <- melt(Zeng_results) %>% tibble()
Zeng_results_melted %<>% rename(marker_set = variable)
Zeng_results_melted %<>% rename(auroc = value)
Zeng_results_melted

Zeng_results_melted %<>% mutate(marker_set = gsub("_markers", "", marker_set))
Zeng_results_melted %<>% mutate(cluster_name = gsub("_cluster", "", cluster_name))

Zeng_results_melted %<>% group_by(marker_set) %>% mutate(auroc_rank = rank(auroc))

max(Zeng_results_melted$auroc_rank) 

self_markers <- Zeng_results_melted %>% filter(cluster_name == marker_set)

self_markers %>% pull(auroc) %>% mean()
self_markers %>% pull(auroc_rank) %>% mean()

#get top hit for each marker set
best_clusters_per_marker <- Zeng_results_melted %>% group_by(marker_set) %>% arrange(-auroc) %>% summarize(best_cluster = first(cluster_name), best_AUROC = first(auroc),
                                                                                                           second_cluster = nth(cluster_name, 2), second_AUROC = nth(auroc,2))
best_clusters_per_marker %<>% mutate(marker_source = if_else(grepl("^[0-9]+$", marker_set), "Zeng", "Macosko"))
best_clusters_per_marker %<>% mutate(best_cluster_source = if_else(grepl("^[0-9]+$", best_cluster), "Zeng", "Macosko"))

best_clusters_per_marker %>% filter(marker_set == "2215")

best_clusters_per_marker %>% filter(best_cluster == marker_set) %>% nrow()
best_clusters_per_marker %>% pull(best_AUROC) %>% mean()


####################################
####################################
####################################
#Self versus best
Zeng_self_results <- read_csv(here("results", "marker_results", "Zeng_self_vs_best.1713323744.csv.gz"))
Zeng_self_results %>% filter(max_AUROC < self_1_vs_all_AUROC) %>% nrow()/ nrow(Zeng_self_results)

Macosko_self_results <- read_csv(here("results", "marker_results", "Macosko_self_vs_best.1713379407.csv.gz"))
Macosko_self_results %>% filter(max_AUROC < self_1_vs_all_AUROC) %>% nrow()/nrow(Macosko_self_results)
Macosko_self_results %>% mutate(source="Macosko")
self_results <- bind_rows(Macosko_self_results %>% mutate(source="Macosko"), Zeng_self_results %>% mutate(source = "Zeng", target_cluster = as.character(target_cluster)))
self_results %>% group_by(source) %>% summarize(n = n(), 
                                                mean_self_vs_best_AUROC = mean(self_vs_best_AUROC),
                                                enriched_self = sum(self_vs_best_AUROC > 0.5),
                                                unique_highest_non_self_cluster = length(unique(highest_non_self_cluster)))

ggplot(data = self_results, aes(x = self_vs_best_AUROC)) + geom_histogram() + theme_bw() + facet_wrap(.~source)


self_results %>% filter(max_AUROC < self_1_vs_all_AUROC) %>% group_by(source) %>% summarize(n = n(), 
                                                                                            mean_self_vs_best_AUROC = mean(self_vs_best_AUROC),
                                                                                            enriched_self = sum(self_vs_best_AUROC > 0.5),
                                                                                            unique_highest_non_self_cluster = length(unique(highest_non_self_cluster)))


for_violin <- bind_rows(self_results %>% select(target_cluster, source, auroc = self_1_vs_all_AUROC) %>% mutate(type = "self_1_vs_all_AUROC"),
                        self_results %>% select(target_cluster, source, auroc = self_vs_best_AUROC) %>% mutate( type = "self_vs_best_AUROC"))

for_violin %<>% mutate(color = recode(source, Macosko = '#E21E25', Zeng = '#4450A2', Siletti = "#E21E25"))
for_violin %<>% mutate(source_label = recode(source, Macosko = 'Nuclei', Zeng = 'Cells'))
for_violin %<>% mutate(Comparison = if_else(type == "self_1_vs_all_AUROC", "Target cluster vs all", "Target vs best\noff-target cluster"))

best_clusters_per_marker %>% filter(best_cluster == marker_set) %>% mutate(marker_set = paste0("Zeng|", marker_set)) %>% filter(marker_set %in% recip_types) %>% 
  nrow()

pdf(here("results", "marker_results","AUROCs_histograms.pdf"), height=6, width = 8)
ggplot(for_violin, aes(x=auroc, fill = Comparison)) + 
  geom_histogram(data= for_violin %>% filter(type == "self_1_vs_all_AUROC")) + 
  geom_histogram(data= for_violin %>% filter(type == "self_vs_best_AUROC")) + 
  facet_wrap(. ~ Comparison) + theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) 
dev.off()


#split by atlas, little difference
ggplot(for_violin, aes(x=auroc, fill = Comparison)) + 
  geom_histogram(data= for_violin %>% filter(type == "self_1_vs_all_AUROC")) + 
  geom_histogram(data= for_violin %>% filter(type == "self_vs_best_AUROC")) + 
  facet_wrap(source_label ~ Comparison) + theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) 
