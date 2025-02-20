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
Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)

base_folder <- here("results/full_run_ZengAWS.1698256033/") #2009 reciprocal result

#CSV files generated in python with examine_region_overlap notebook
nuclei_region_counts <- read_csv(here("results", "region_profile_compare", "Macosko_region_counts.csv"))
nuclei_region_counts %<>% filter(cell_count_Macosko != 0)
nuclei_region_counts %<>% group_by(region_remap) %>% summarize(cluster_count = n(), cell_count = sum(cell_count_Macosko), median_cells_per_cluster = median(cell_count_Macosko))
nuclei_region_counts %<>% mutate(Dataset = "Nuclei")

cell_region_counts <- read_csv(here("results", "region_profile_compare", "Zeng_region_counts.csv"))
cell_region_counts %<>% filter(cell_count_Zeng != 0)
cell_region_counts %<>% group_by(region_remap) %>% summarize(cluster_count = n(), cell_count = sum(cell_count_Zeng), median_cells_per_cluster = median(cell_count_Zeng))
cell_region_counts %<>% mutate(Dataset = "Cells")

joined <- inner_join(cell_region_counts, nuclei_region_counts, by="region_remap", suffix = c("_cells", "_nucleus"))
joined %>% select(region_remap, cluster_count_cells, cluster_count_nucleus)
joined %>% select(region_remap, cell_count_cells, cell_count_nucleus)
joined %>% write_csv(here("results", "region_profile_compare", "joined_region_counts.csv"))

colnames(joined)
region_order <- joined %>% select(region_remap, cell_count_nucleus) %>% 
  arrange(-cell_count_nucleus) %>% pull(region_remap) %>% as.character() %>% unique()
joined %<>% mutate(region_remap = factor(region_remap, levels = region_order))

#colors and labels were changed 
ggplot(joined, aes(x = region_remap, y = cell_count_nucleus, fill = region_remap)) + 
  geom_bar(stat="identity") + ylim(-7e5, 7e5) + theme_bw() +
  theme(legend.position = 'none')

ggplot(joined, aes(x = region_remap, y = cell_count_cells, fill = region_remap)) + 
  geom_bar(stat="identity") + ylim(-7e5, 7e5) + theme_bw() +
  theme(legend.position = 'none')

#cluster counts
ggplot(joined, aes(x = region_remap, y = cluster_count_nucleus, fill = region_remap)) + 
  geom_bar(stat="identity") + ylim(-3000, 3000) + theme_bw() +
  theme(legend.position = 'none')

ggplot(joined, aes(x = region_remap, y = cluster_count_cells, fill = region_remap)) + 
  geom_bar(stat="identity") + ylim(-3000, 3000) + theme_bw() +
  theme(legend.position = 'none')


#Cumulative plot
#uses output from a full python metaneighbor run via run_on_full_merged_ZengAWS (base_folder variable above)
#saved output is provided in full_run_ZengAWS.1698256033
merged_obs <- read_csv(paste0( base_folder, "merged.obs.csv.zip"))

paste0("Number of genes used:", read_csv(paste0( base_folder, "merged.var.csv.zip")) %>% nrow())

#cell counts per study
merged_obs %>% group_by(study_id) %>% count()

merged_obs %<>% rename(cell_id = "...1")
cell_type_sizes <- merged_obs %>% group_by(cell.type, study_id) %>% count()
cell_type_sizes %<>% ungroup()
cell_type_sizes %<>% mutate(cell.type = paste0(study_id, "|", cell.type))

for_plot <- cell_type_sizes %>% group_by(study_id) %>% arrange(-n) %>% 
  mutate(cumulative_cell_fraction = cumsum(n)/sum(n), cumulative_cluster_fraction  = row_number()/n())
for_plot %<>% bind_rows(tibble(study_id = "Zeng", cumulative_cell_fraction = 0, cumulative_cluster_fraction =0))
for_plot %<>% bind_rows(tibble(study_id = "Macosko", cumulative_cell_fraction = 0, cumulative_cluster_fraction =0))
for_plot %<>% mutate(Dataset = if_else(study_id == "Zeng", "Cells", "Nuclei"))

#pdf(paste0(base_folder, "cluster_size_gini_plot.pdf"), height=5, width = 5)
ggplot(data = for_plot, aes(y=cumulative_cell_fraction, x=cumulative_cluster_fraction, color = Dataset)) + 
  geom_line() + theme_bw() + coord_fixed(expand=FALSE) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + 
  ylab("Cumulative cell fraction") + xlab("Cumulative cluster fraction") +
  theme(legend.position = "bottom") + scale_color_manual(values = c("#4450A2", "#E21E25"))
#dev.off()  

