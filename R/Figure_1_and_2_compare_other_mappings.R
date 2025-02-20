library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(readr)
library(here)
library(magrittr)
library(conflicted)
library(eulerr)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)

#for writing files out
base_folder <- here("results","full_run_ZengAWS.1698256033/") #updated Allen data - 2009 recips

Zeng_cluster_info_new <- read_tsv(here("data", "whole_mouse_brain", "zeng", "from_aws", "AIT21.0", "AIT21_annotation_freeze_081523.tsv"))
Zeng_cluster_info_new %<>% mutate(cl = as.character(cl))
colnames(Zeng_cluster_info_new)
LQ_cluster_cl <- Zeng_cluster_info_new %>% filter(subclass_label == "LQ") %>% pull(cl)
LQ_cluster_ids <- Zeng_cluster_info_new %>% filter(subclass_label == "LQ") %>% pull(cluster_id)
Zeng_cluster_info_new %<>% filter(subclass_label != "LQ")
all_Zeng_clusters <- unique(Zeng_cluster_info_new %>% pull(cl))
all_Zeng_cluster_id <- unique(Zeng_cluster_info_new %>% pull(cluster_id))

length(all_Zeng_clusters)

Liu_matches <- read_csv(here("data", "Hanqing_Liu_etal", "41586_2023_6805_MOESM7_ESM.csv"))
#pretty sure these are cl from an older version - WB.annotation.cl.Level2.20221011.csv - in that microglia is the same ID
#we need to filter it for CL's that remain
nrow(Liu_matches)
Liu_matches %<>% select(mC_Cluster = `mC Cluster`, matched_cluster = `Matched AIBS 10X RNA Clusters` ) 
Liu_matches %>% filter(matched_cluster == "5283")
Liu_matches %>% filter(matched_cluster == "")
Liu_matches %<>% separate_rows(matched_cluster, sep = ", ") 
#lets remove the Liu matches that are no longer in the Zeng annotations - the cl values seem stable
nrow(Liu_matches)
Liu_matches %>% filter(!matched_cluster %in% all_Zeng_clusters) %>% pull(matched_cluster) %>% unique()
Liu_matches_universe_filtered <- Liu_matches %>% filter(matched_cluster %in% all_Zeng_clusters)
nrow(Liu_matches_universe_filtered)

Liu_matches %>% group_by(mC_Cluster) %>% count() %>% group_by(n) %>% count()

#just use the ones that are in the universe
Liu_match_ids <- Liu_matches_universe_filtered %>% pull(matched_cluster) %>% unique()
length(Liu_match_ids)

#how many 1-1 mappings are there?
mapped_to_one_Zeng <- Liu_matches %>% group_by(mC_Cluster) %>% count() %>% filter(n == 1) %>% pull(mC_Cluster)
length(mapped_to_one_Zeng)
mapped_to_one_Liu <- Liu_matches %>% group_by(matched_cluster) %>% count() %>% filter(n == 1) %>% pull(matched_cluster)
length(mapped_to_one_Liu)
Liu_matches %>% filter(mC_Cluster %in%  mapped_to_one_Zeng, matched_cluster %in% mapped_to_one_Liu) %>% nrow() #145
Liu_matches %>% pull(mC_Cluster) %>% unique() %>% length()


#only male mice cells used
Songpeng_matches <- read.csv(here("data", "Songpeng_Zu_etal", "SI Table 5 Transfer label scores for the integration of the snATAC-seq with the scRNA-seq data. Both cluster and subclass level scores are presented", "sa2.all.cl2L4.transferLabelScore.csv"))
dim(Songpeng_matches)
Songpeng_matches$Zeng_id = rownames(Songpeng_matches)
Songpeng_matches_melted <- melt(Songpeng_matches, id.vars = "Zeng_id") %>% tibble()
Songpeng_matches_melted %>% group_by(variable) %>% summarize(n=n(), sum_values = sum(value)) #check that the rows add up to one - proportion of cells matched

#from Songpeng et al:
#"For each L4-level subtype, we used the corresponding top 3 clusters in the scRNA-seq data as the 
#candidate annotations, then mapped the three clusters to the subclasses defined in the scRNA-seq data, 
#and manually checked whether they were consistent on mouse brain major regions and gene markers."
#to get 1040 ID's in Extended data Fig 6 g+h
Songpeng_matches_melted %<>% group_by(variable) %>% arrange(-value) %>% filter(row_number() < 2)
Songpeng_matches_melted %>% arrange(-value) %>% filter(value > 0.5)
cells_matched_in_top_one <- Songpeng_matches_melted %>% group_by(variable) %>% summarize(n=n(), sum_values = sum(value))
cells_matched_in_top_one %>% pull(sum_values) %>% median() #proportion matched to the top three
cells_matched_in_top_one %>% pull(sum_values) %>% mean()
Songpeng_matches_melted %>% pull(Zeng_id) %>% unique() %>% length()
Songpeng_matches_melted %>% pull(Zeng_id) %>% unique() %>% length()/5322

Songpeng_candidate_ids <- Songpeng_matches_melted %>% pull(Zeng_id) %>% unique()
length(Songpeng_candidate_ids)


#related image: https://www.nature.com/articles/s41586-023-06817-8/figures/14
Carla_Winter_table <- read_xlsx(here("data", "Carla_Winter_etal","Supp Table_5_He_AIBS_Mapping_Cluster_Summary.xlsx")) %>% 
  mutate(AIBS_cluster_id = gsub(" .*", "", AIBS_cluster_label)) %>%
  mutate(AIBS_cluster_id = gsub("^0*", "", AIBS_cluster_id)) 
Carla_Winter_table %>% pull(He_label) %>% unique() %>% length()
#how many 1-1
mapped_to_one_Zeng <- Carla_Winter_table %>% group_by(He_label) %>% count() %>% filter(n == 1) %>% pull(He_label) #He is the last author
length(mapped_to_one_Zeng)
mapped_to_one_He <- Carla_Winter_table %>% group_by(AIBS_cluster_id) %>% count() %>% filter(n == 1) %>% pull(AIBS_cluster_id)
length(mapped_to_one_He)
Carla_Winter_table %>% filter(He_label %in% mapped_to_one_Zeng, AIBS_cluster_id %in% mapped_to_one_He) #zero

#convert to cl ID's from cluster_id's
Carla_Winter_table %<>% left_join(Zeng_cluster_info_new %>% select(cl, AIBS_cluster_label = cluster_id_label))

Carla_Winter_ids <- Carla_Winter_table %>% pull(cl) %>% unique()
length(Carla_Winter_ids)


# Create a list of sets
sets <- list(`Winter et al.` = Carla_Winter_ids,
             `Songpeng et al.` = Songpeng_candidate_ids,
             `Liu et al.` = Liu_match_ids, 
             all=all_Zeng_clusters)

# Create a Venn diagram
fit <- euler(sets)

# Plot the Venn diagram
pdf(here("results", "Euler_diagram_comparison_other_mappings.pdf"), height=6, width = 6)
plot(fit, quantities = TRUE)
dev.off()
plot(fit)

length(Liu_match_ids)
length(Songpeng_candidate_ids)
intersect(Liu_match_ids, Songpeng_candidate_ids) %>% length() #p-value computed with https://brain.shinyapps.io/hyper/
#all three mapped
intersect(intersect(Carla_Winter_ids, Songpeng_candidate_ids), Liu_match_ids) %>% length()

#calculate a p-value on getting the intersection of all three
iterations = 1000
all_types <- 1:5322
all_resuts <- tibble()
for(iteration in 1:iterations) {
  a <- sample(all_types, length(Carla_Winter_ids))
  b <- sample(all_types, length(Songpeng_candidate_ids))
  c <- sample(all_types, length(Liu_match_ids))
  
  all_resuts %<>% bind_rows(tibble(iteration, intersect = intersect(intersect(a, b), c) %>% length()))
}
#significantly depleted for matches using all
all_resuts %>% pull(intersect) %>% mean()
all_resuts %>% filter(intersect >= 32) %>% nrow() / iterations
all_resuts %>% filter(intersect <= 32) %>% nrow() / iterations

#intersect the reciprocal hits for the second figure
top_hits <- read_csv(paste0(base_folder, "top_hits.0.95.csv"))

recip_types <- union(top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_1`),
                     top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_2`))

#needs a list of all types for this
if(file.exists(paste0(base_folder, "merged.obs.csv"))) {
  merged_obs <- read_csv(paste0( base_folder, "merged.obs.csv"))
} else {
  merged_obs <- read_csv(paste0( base_folder, "merged.obs.csv.zip"))
}
merged_obs %<>% rename(cell_id = "...1") 
merged_obs %<>% mutate(cell.type = paste0(study_id, "|", cell.type))

recip_overlaps <- merged_obs %>% select(cell.type, study_id) %>% mutate() %>% distinct() %>% 
  mutate(is_recip = cell.type %in% recip_types) %>% group_by(study_id, is_recip) %>% count() 
recip_overlaps %<>% mutate(Dataset = if_else(study_id=="Macosko", "Nucleus", "Cell")) %>% ungroup() %>% select(-study_id)

slide_seq_mapped <- read_delim(here("data", "whole_mouse_brain", "macosko", "braincelldata.org", "Mapping_matrices", "Puck_Num_01.mapping.MappedCellTypes.txt"), delim="=", col_names = F)  %>% pull(X2)
recip_types_short <- gsub("Zeng.", "", recip_types)
recip_types_short %<>% gsub("Macosko.", "", .)
recip_overlaps %<>% bind_rows(tibble(Dataset = "Slideseq\nmapped", is_recip = FALSE, n = length(setdiff(slide_seq_mapped, recip_types_short))))
recip_overlaps %<>% bind_rows(tibble(Dataset = "Slideseq\nmapped", is_recip = TRUE, n = length(base::intersect(slide_seq_mapped, recip_types_short))))

#pvalues from https://brain.shinyapps.io/hyper/
#X >= 953, p =  2.354e-26 
recip_overlaps %<>% bind_rows(tibble(Dataset = "Liu et al.", is_recip = FALSE, n = length(setdiff(Liu_match_ids, recip_types_short))))
recip_overlaps %<>% bind_rows(tibble(Dataset = "Liu et al.", is_recip = TRUE, n = length(base::intersect(Liu_match_ids, recip_types_short))))
#X <= 1639, p =  9.141e-26  #worse than expected

recip_overlaps %<>% bind_rows(tibble(Dataset = "Songpeng et al.", is_recip = FALSE, n = length(setdiff(Songpeng_candidate_ids, recip_types_short))))
recip_overlaps %<>% bind_rows(tibble(Dataset = "Songpeng et al.", is_recip = TRUE, n = length(base::intersect(Songpeng_candidate_ids, recip_types_short))))

recip_overlaps %<>% bind_rows(tibble(Dataset = "Winter et al.", is_recip = FALSE, n = length(setdiff(Carla_Winter_ids, recip_types_short))))
recip_overlaps %<>% bind_rows(tibble(Dataset = "Winter et al.", is_recip = TRUE, n = length(base::intersect(Carla_Winter_ids, recip_types_short))))

##of the 32 mapped in all three - 17 are recips
base::intersect(Liu_match_ids, base::intersect(Carla_Winter_ids, Songpeng_candidate_ids))
base::intersect(Liu_match_ids, base::intersect(Carla_Winter_ids, Songpeng_candidate_ids)) %>% length()
base::intersect(recip_types_short, base::intersect(Liu_match_ids, base::intersect(Carla_Winter_ids, Songpeng_candidate_ids))) %>% length()
#X >= 48, p =  0.002 
length(Carla_Winter_ids)

recip_overlaps %<>% mutate(Dataset = factor(Dataset, levels = c("Cell", "Nucleus", "Slideseq\nmapped", "Liu et al.", "Songpeng et al.", "Winter et al.")))


pdf(paste0(base_folder, "recip_overlaps.", top_hits %>% nrow(), "_pairs.pdf"), height=5, width = 7)
recip_overlaps %<>% mutate(Clusters = if_else(is_recip, "Reciprocal hit", "All"))
ggplot(data = recip_overlaps, aes(x = Dataset, y = n, fill = Clusters)) + geom_bar(stat = "identity") + 
  theme_bw() + xlab("") + scale_fill_manual(values=c("grey", "black")) + ylab("Clusters")
dev.off()

#calculate proportions to add to the figure and use in the text
recip_overlaps %<>% group_by(Dataset) %>% mutate(dataset_n = sum(n), prop = n/dataset_n)
recip_overlaps
