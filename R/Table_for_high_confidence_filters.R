library(ggplot2)
library(tidyr)
library(conflicted)
library(magrittr)
library(reshape2)
library(dplyr)
library(here)
library(readr)
conflict_prefer("count", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)

base_folder <- here("results/full_run_ZengAWS.1698256033/") #2009 hits

top_hits <- read_csv(paste0(base_folder, "top_hits_asymmetric.0.99.csv"))
top_hits %<>% filter(Match_type == "Reciprocal_top_hit")

split_references <- bind_rows(top_hits %>% select(reference = `Study_ID|Celltype_1`, paired_reciprocal = `Study_ID|Celltype_2`),
                              top_hits %>% select(reference = `Study_ID|Celltype_2`, paired_reciprocal = `Study_ID|Celltype_1`))

split_references %<>% mutate(ref_study = getStudyId(reference))

#check if the paired_reciprocal is in the best vs next
#get the AUROC for that
one_vs_one_AUROCs <- read_csv(paste0(base_folder, "aurocs_1v1.csv.gz"), guess_max = 200000000)
one_vs_one_AUROCs %<>% rename(target = "...1")
#melting to pull out more details
one_vs_one_AUROCs_melted <- melt(one_vs_one_AUROCs, na.rm = T, value.name= "auroc", variable.name = "reference") %>% tibble()
one_vs_one_AUROCs_melted %<>% mutate(reference = as.character(reference), target = as.character(target))
one_vs_one_AUROCs_melted %<>% mutate(reference_study = getStudyId(reference), target_study = getStudyId(target))
#add in sizes of targets
one_vs_one_AUROCs_melted 
one_vs_one_AUROCs_melted %>% group_by(reference) %>% count()

split_references %<>% left_join(one_vs_one_AUROCs_melted %>% select(reference, paired_reciprocal = target, best_vs_next_auroc = auroc))

split_references %>% group_by(is.na(best_vs_next_auroc)) %>% count() #45
split_references %>% filter(is.na(best_vs_next_auroc))

#remove those that didn't show up in the best vs next matchups
split_references %<>% filter(!is.na(best_vs_next_auroc))  

#how many pairs do we have at this point?
top_hits %>% filter(`Study_ID|Celltype_1` %in% (split_references %>% pull(reference)),
                    `Study_ID|Celltype_2` %in% (split_references %>% pull(reference))) %>% nrow()

#filter for both ends of surviving 
split_references %>% nrow()
split_references %<>% filter(best_vs_next_auroc > 0.6)
split_references %>% nrow()
hist(split_references %>% pull(best_vs_next_auroc))
split_references %>% pull(best_vs_next_auroc) %>% median()

top_hits %<>% filter(`Study_ID|Celltype_1` %in% (split_references %>% pull(reference)),
                     `Study_ID|Celltype_2` %in% (split_references %>% pull(reference)))
top_hits %>% nrow()

#order so Zeng is first
top_hits %<>% rowwise() %>% mutate(first = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[2],
                                   second = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[1])
top_hits %<>% mutate(`Study_ID|Celltype_1` = first, `Study_ID|Celltype_2` = second)
top_hits %<>% select(-first, -second)

top_hits %>% write_csv(paste0(base_folder, "top_hits_asymmetric.0.99.best_vs_next.0.6.filtered.csv"))

#compare to the 2009
top_hits_symmetric <- read_csv(paste0(base_folder, "top_hits.0.95.csv"))
top_hits_symmetric %<>% filter(Match_type == "Reciprocal_top_hit")

top_hits_symmetric %<>% rowwise() %>% mutate(first = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[2],
                                             second = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[1])
top_hits_symmetric %<>% mutate(`Study_ID|Celltype_1` = first, `Study_ID|Celltype_2` = second)
top_hits_symmetric %<>% select(-first, -second)


inner_join(top_hits_symmetric %>% select(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`),
           top_hits %>% select(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`)) %>% nrow()

anti_join(top_hits %>% select(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`),
          top_hits_symmetric %>% select(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`)) %>% as.data.frame()
