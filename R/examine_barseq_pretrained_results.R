library(readxl)
library(tidyr)
library(MetaNeighbor)
library(reshape2)
library(magrittr)
library(dplyr)
library(here)
library(readr)
library(tibble)
library(conflicted)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)

base_barseq_folder <- here("data/whole_mouse_brain/processed/barseq/")                           

base_folder <- here("results/full_run_ZengAWS.1698256033/") #2009 hits location

pretrained_results <- read_csv(here("results/barseq_109_pretrained_R.1728586264/Zeng_pretrained_aurocs.csv.gz") )
dim(pretrained_results)

n_BARseq <- nrow(pretrained_results)
n_Zeng <- ncol(pretrained_results)

pretrained_results %<>% bind_rows(read_csv(here("results/barseq_109_pretrained_R.1728586264/Macosko_pretrained_aurocs.csv.gz")))
n_Macosko <- ncol(pretrained_results) - n_Zeng

pretrained_results %<>% melt() %>% tibble()

pretrained_results %<>% rename(target = `...1`)
pretrained_results %<>% mutate(study = getStudyId(as.character(variable)))

#get the top hit per study and type
pretrained_results %<>% group_by(target, study) %>% slice_max(value)

#test for top hits
top_hits <- read_csv(paste0(base_folder, "top_hits.0.95.csv"))

top_hits %<>% filter(Match_type == "Reciprocal_top_hit")
top_hits %>% nrow()
recip_clusters <- c(top_hits %>% pull(`Study_ID|Celltype_1`), top_hits %>% pull(`Study_ID|Celltype_2`)) %>% unique()
length(recip_clusters)

pretrained_results %<>% mutate(is_recip_top_hit = variable %in% recip_clusters)
pretrained_results %>% pull(variable) %>% unique() %>% length()
pretrained_results %>% filter(is_recip_top_hit)
pretrained_results %>% group_by(is_recip_top_hit) %>% count()
pretrained_results_before_summary <- pretrained_results

pretrained_results %>% group_by(study) %>% summarize(mean_AUROC = mean(value))
pretrained_results %>% group_by(study) %>% summarize(median_AUROC = median(value))

pretrained_results %<>% group_by(is_recip_top_hit, study) %>% count() %>% mutate(proportion = n/n_BARseq) %>% arrange(study)
pretrained_results$n_original = c(n_Macosko-nrow(top_hits),nrow(top_hits),n_Zeng-nrow(top_hits),nrow(top_hits))

pretrained_results %<>% pivot_wider(names_from = is_recip_top_hit, values_from = c(n, proportion, n_original))
pretrained_results %<>% mutate(proportion_TRUE_original = n_original_TRUE/(n_original_TRUE+n_original_FALSE))

#calculate odds ratios
pretrained_results %>% select(study, proportion_TRUE, proportion_TRUE_original) %>%
  mutate(
    odds_ratio = (proportion_TRUE / (1 - proportion_TRUE)) / 
      (proportion_TRUE_original / (1 - proportion_TRUE_original))
  )

#calculate p-values
pretrained_results %>%
  rowwise() %>%
  mutate(
    odds_ratio = (proportion_TRUE / (1 - proportion_TRUE)) / 
      (proportion_TRUE_original / (1 - proportion_TRUE_original)),
    p_value = fisher.test(matrix(c(n_TRUE, n_FALSE, n_original_TRUE, n_original_FALSE), nrow = 2))$p.value
  )

pretrained_results <- pretrained_results_before_summary
#97 - clusters with any recip hit
pretrained_results %>% group_by(target) %>% filter(is_recip_top_hit) %>% count() %>% filter(n>0) %>% nrow()

#need to find the number of clusters that map to both ends of a reciprocal pair then compare to curated file
targets_with_both_recip <- pretrained_results %>% group_by(target) %>% filter(is_recip_top_hit) %>% count() %>% filter(n==2) %>% pull(target)
length(targets_with_both_recip)
#how many of those pairs are recip pairs in the tophits? need that filter too
pretrained_results_recip_pairs <- pretrained_results %>% group_by(target) %>% filter(is_recip_top_hit) %>% 
  arrange(study) %>% summarize(n=n(), Mac_hit = first(variable), Zeng_hit = last(variable)) %>% filter(n==2) 

#those that don't pass
anti_join(pretrained_results_recip_pairs, top_hits %>% select(Zeng_hit = `Study_ID|Celltype_1`, Mac_hit = `Study_ID|Celltype_2`))
pretrained_results_recip_pairs %<>% inner_join(top_hits %>% select(Zeng_hit = `Study_ID|Celltype_1`, Mac_hit = `Study_ID|Celltype_2`))
#key number - both sides of a recip pair
nrow(pretrained_results_recip_pairs)

#check that those used in the figure also there
pretrained_results_recip_pairs %>% filter(target == "BARseq|SubCTX_4")
pretrained_results_recip_pairs %>% filter(target == "BARseq|SubCTX_7")
pretrained_results_recip_pairs %>% filter(target == "BARseq|P D?")
pretrained_results_recip_pairs %>% pull(target) %>% sort()

#write out curation spreadsheet/suppelment table for annotation 
region_check_table <- pretrained_results_recip_pairs
region_check_table %<>% rename(Zeng_name = Zeng_hit, Macosko_name = Mac_hit)
region_check_table %<>% mutate(Macosko_name = getCellType(Macosko_name), Zeng_name = getCellType(Zeng_name))
region_check_table %<>% select(-n)
region_check_table %<>% rename(name = target)
top_region_barseq <- read_csv(paste0(base_barseq_folder, "top_region_per_cluster.csv"))

ccfv3_acros <- read_xlsx(here("data" ,"Allen_CCFv3","1-s2.0-S0092867420304025-mmc2.xlsx"), skip = 1)
ccfv3_acros %<>% select(abbreviation, region_name =`full structure name`) %>% distinct()

region_check_table %<>% mutate(name = getCellType(name))
region_check_table %<>% rename(cluster_H3 = name)

#add in top region per barseq cluster
region_check_table %<>% inner_join(top_region_barseq %>% select(cluster_H3 = cluster, barseq_region_top, -n))
region_check_table

#add in dissection region and other data from Macosko
Macosko_updated_annots <- read_tsv(here("data" , "whole_mouse_brain",  "macosko", "from_terra_non_uniform", "CellType_Metadata.tsv"))
colnames(Macosko_updated_annots)
Macosko_updated_annots %<>% select(Macosko_name = Annotation, 
                                   Macosko_top_dissectate_postQC = top_dissectate_postQC, 
                                   Macosko_Imputed_Top_DeepCCF = Imputed_Top_DeepCCF, 
                                   Macosko_Imputed_Top_TopStruct = Imputed_Top_TopStruct, 
                                   Macosko_Imputed_Top_DeepCCF = Imputed_Top_DeepCCF)
Macosko_updated_annots
region_check_table
region_check_table %<>% inner_join(Macosko_updated_annots)
region_check_table

#load in the Macosko ones that made it to slide seq so we can inspect them
mapped_types <- read_delim(here("data","whole_mouse_brain", "macosko", "braincelldata.org", "Mapping_matrices", "Puck_Num_01.mapping.MappedCellTypes.txt"), delim="=", col_names = F) %>% pull(X2)
length(mapped_types)
region_check_table %<>% mutate(is_slideseq_mapped = Macosko_name %in% mapped_types)
region_check_table

Zeng_cluster_info_new <- read_tsv(here("data/whole_mouse_brain/zeng/from_aws/AIT21.0/AIT21_annotation_freeze_081523.tsv"))
Zeng_cluster_info_new %<>% mutate(cl = as.character(cl))
colnames(Zeng_cluster_info_new)
Zeng_cluster_info_new[1,] %>% as.data.frame()
Zeng_cluster_info_new %<>% select(Zeng_name = cl, Zeng_MERFISH_ID = cluster_id, Zeng_super_label = supertype_label, Zeng_anatomical_annotation = anatomical_annotation, Zeng_CCF_broad_freq = CCF_broad.freq)

region_check_table %<>% inner_join(Zeng_cluster_info_new)
region_check_table %<>% rename(BARseq_cluster = cluster_H3)

#this is the curated spreadsheet with all the information added in for the pairs
dir.create(here("results", "barseq_mapping"))
region_check_table %>% write_csv(here("results", "barseq_mapping", "pretrained_precuration_annotated.csv"))

