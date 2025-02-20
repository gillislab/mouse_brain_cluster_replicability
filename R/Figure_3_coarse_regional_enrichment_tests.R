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
conflict_prefer("count", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
here()

Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)

base_folder <- here("results/full_run_ZengAWS.1698256033/") #2009 hits

nuclei_region_counts <- read_csv(here("results", "region_profile_compare", "Macosko_region_counts.csv"))
nuclei_region_counts %<>% select(region_remap, cell_type = ClusterNm, cell_count = cell_count_Macosko)
nuclei_region_counts %<>% mutate(Dataset = "Macosko")

cell_region_counts <- read_csv(here("results", "region_profile_compare", "Zeng_region_counts.csv"))
cell_region_counts %<>% select(region_remap, cell_type = cl, cell_count = cell_count_Zeng)
cell_region_counts %<>% mutate(cell_type = as.character(cell_type))
cell_region_counts %<>% mutate(Dataset = "Zeng")

joined_region_counts <- bind_rows(nuclei_region_counts, cell_region_counts)

top_hits <- read_csv(paste0(base_folder, "top_hits.0.95.csv"))
all_recip_hits <- top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% select(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`) %>% pivot_longer(everything()) %>% pull(value)

joined_region_counts

joined_region_counts %<>% mutate(is_recip_hit = paste0(Dataset, "|", cell_type) %in% all_recip_hits)
joined_region_counts %>% group_by(is_recip_hit, Dataset) %>% summarize(cell_count = sum(cell_count))

joined_region_counts %<>% group_by(region_remap, Dataset, is_recip_hit) %>% summarize(cell_count = sum(cell_count))
joined_region_counts %<>% bind_rows(joined_region_counts %>% group_by(region_remap, is_recip_hit) %>% summarize(cell_count = sum(cell_count)) %>% mutate(Dataset="Combined"))
tail(joined_region_counts)

#table with region, total cells, recip cells and non-recip cells, all recip cells, and study
joined_region_counts %<>% pivot_wider(names_from = is_recip_hit, values_from = cell_count, names_prefix = "recip_hit_")
joined_region_counts %<>% mutate(total_cells_in_region = recip_hit_TRUE + recip_hit_FALSE)
joined_region_counts %<>% mutate(percent_recip = recip_hit_TRUE / total_cells_in_region)
joined_region_counts %>% ungroup() %>% arrange(-percent_recip)
joined_region_counts %>% ungroup() %>% arrange(percent_recip)
#add in totals
joined_region_counts %<>% inner_join(joined_region_counts %>% group_by(Dataset) %>% summarize(dataset_recip_hit_TRUE = sum(recip_hit_TRUE),dataset_recip_hit_FALSE = sum(recip_hit_FALSE), dataset_total_cells = sum(total_cells_in_region)))

joined_region_counts %<>% mutate(dataset_prop_recip = dataset_recip_hit_TRUE/ dataset_total_cells)

# then compute - one sided p-value twice?
joined_region_counts %<>% mutate(p_enriched = phyper(recip_hit_TRUE, dataset_recip_hit_TRUE, dataset_recip_hit_FALSE, total_cells_in_region, lower.tail = FALSE)  + 
                                   dhyper(recip_hit_TRUE, dataset_recip_hit_TRUE, dataset_recip_hit_FALSE, total_cells_in_region))

joined_region_counts %<>% mutate(p_depleted = phyper(recip_hit_TRUE, dataset_recip_hit_TRUE, dataset_recip_hit_FALSE, total_cells_in_region, lower.tail = TRUE))
joined_region_counts %<>% mutate(difference_from_dataset_average = percent_recip - dataset_prop_recip)
joined_region_counts %>% write_csv(here(base_folder, "top_hit_enrichments", "regional_enrichment_results.csv"))

joined_region_counts

joined_region_counts %<>% mutate(Dataset= if_else(Dataset=="Zeng", "Cells", Dataset))
joined_region_counts %<>% mutate(Dataset= if_else(Dataset=="Macosko", "Nuclei", Dataset))

joined_region_counts %<>% ungroup()

#add in MERFISH data from examine_MERFISH_calls_for_recips.R
by_region <- read_csv(here(base_folder, "top_hit_enrichments", paste0("MERFISH_regional_enrichment_results.", top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% nrow(), "_pairs.csv")))
joined_region_counts %<>% bind_rows(by_region %>% rename(region_remap = parcellation_division_remap))
joined_region_counts %<>% mutate(Dataset = factor(Dataset, levels = rev(c("Cells", "Nuclei", "Combined", "MERFISH"))))

pdf(here(base_folder, "top_hit_enrichments", "dissection_region_heatmap_plus_MERFISH_no_combined.pdf"), height=3.5, width = 8)
max_diff <- joined_region_counts %>% pull(difference_from_dataset_average) %>% abs() %>% max()
color_mapping <- rev((grDevices::colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu")))(100))

new_ordering <- joined_region_counts %>% group_by(region_remap) %>% summarize(difference_from_dataset_average = mean(difference_from_dataset_average))
new_ordering %<>% arrange(-difference_from_dataset_average)

#max diff is unchanged, sorted by average
joined_region_counts %<>% mutate(region_remap = factor(region_remap, levels = new_ordering %>% pull(region_remap)))
joined_region_counts %<>% filter(Dataset != "Combined")
ggplot(data=joined_region_counts %>% filter(!is.na(region_remap)), aes(x=region_remap, y=Dataset, fill = difference_from_dataset_average) ) +
  geom_tile() + scale_fill_gradientn(colors = color_mapping, limits = c(-max_diff, max_diff)) + theme_bw() +
  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) + ylab("") + labs(fill = "Change in proportion of\ncells in reciprocal clusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

