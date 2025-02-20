library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(here)
library(magrittr)
library(conflicted)
conflict_prefer("first", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
conflicts_prefer(base::union)

base_folder <- here("results/full_run_ZengAWS.1698256033/") #updated Allen data - 2009 recips
h5ad_base <- read_csv(here("results", "MERFISH500", "H5ad_calls.csv.gz")) #written out by centroid_correlations.bigcat.R
h5ad_base %<>% mutate(slice = as.numeric(gsub(".*[.]", "", brain_section_label)))
h5ad_base

#uncomment/comment the below lines to generate data for the different top hit sets
top_hits <- read_csv(paste0(base_folder, "top_hits.0.95.csv"))

#top_hits <- read_csv(paste0(base_folder, "top_hits_asymmetric.0.99.best_vs_next.0.6.filtered.csv"))

#Zeng_markers_on_Macosko_one_to_one_mapping.csv
#top_hits <- read_csv(here("results/marker_results/Top_hit_format_for_Zeng_markers_on_Macosko_one_to_one_mapping.csv"))


if ("...1" %in% colnames(top_hits)) { top_hits %<>% select(-`...1`) }
top_hits %<>% rowwise() %>% mutate(first = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[2],
                                   second = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[1])
top_hits %<>% mutate(`Study_ID|Celltype_1` = first, `Study_ID|Celltype_2` = second)
top_hits %<>% select(-first, -second)
top_hits %>% group_by(Match_type) %>% count()
top_hits %<>% filter(Match_type == "Reciprocal_top_hit")
top_hits %>% nrow()
top_hits %<>% mutate(merged_name = paste0(`Study_ID|Celltype_1`, "__",  `Study_ID|Celltype_2`))

Zeng_calls_direct <- read_csv(here("data","whole_mouse_brain",  "zeng", "MERFISH-C57BL6J-638850", "20230830", "cell_metadata.csv"))
Zeng_calls_direct %>% pull(cluster_alias) %>% max() #confirm it's in the cl ID space
Zeng_calls_direct %<>% mutate(cluster_alias = as.character(cluster_alias))
Zeng_calls_direct %>% pull(cluster_alias) %>% unique() %>% length()
Zeng_calls_direct %>% pull(average_correlation_score) %>% min() #confirm the threshold
Zeng_calls_direct %<>% mutate(cluster_alias = paste0("Zeng|", cluster_alias))
Zeng_calls_direct %<>% left_join(top_hits %>% select(cluster_alias = `Study_ID|Celltype_1`, merged_name))

#add in merged name
Zeng_calls_direct %<>% select(cell_label, x, y, z, Celltype_Zeng = cluster_alias, Zeng_merged_name = merged_name)
#join with h5ad to cover all cells
joined <- left_join(h5ad_base, Zeng_calls_direct)

Macosko_calls_direct <- read_csv(here("results/MERFISH500/Macosko_predictions/joined_scores.df.full.csv"))
Macosko_calls_direct %>% pull(prob) %>% min()
Macosko_calls_direct %<>% filter(average_correlation_score > 0.5) #mirroring the Zeng threshold
Macosko_calls_direct %>% pull(prob) %>% min()
Macosko_calls_direct %>% pull(prob) %>% hist()
Macosko_calls_direct %>% pull(average_correlation_score) %>% min()
Macosko_calls_direct %<>% rename(cell_label = cell_id)

Macosko_calls_direct %<>%
  separate(pred.cl, into = c("id", "Celltype_Macosko"), sep = "=")

Macosko_calls_direct %<>% mutate(Celltype_Macosko = paste0("Macosko|", Celltype_Macosko))
Macosko_calls_direct %<>% left_join(top_hits %>% select(Celltype_Macosko = `Study_ID|Celltype_2`, merged_name))
Macosko_calls_direct %>% filter(!is.na(merged_name)) 

Macosko_calls_direct %<>% select(cell_label, Celltype_Macosko, Macosko_merged_name = merged_name)

joined <- left_join(joined, Macosko_calls_direct)
joined %>% group_by(Zeng_merged_name == Macosko_merged_name) %>% count()
joined %>% pull(Celltype_Macosko) %>% unique() %>% length()
joined %>% pull(Celltype_Zeng) %>% unique() %>% length()
joined %>% filter(!is.na(Celltype_Zeng)) %>% nrow()
joined %>% filter(!is.na(Celltype_Macosko)) %>% nrow()

#add in CCF annotations from the API
ccf_annotations <- read_csv(here("data", "whole_mouse_brain", "zeng", "from_API", "ccf_coordinates_MERFISH-C57BL6J-638850.csv.gz"))
base::intersect(joined$cell_label, ccf_annotations$cell_label) %>% length()
ccf_annotations %>% group_by(parcellation_category) %>% count() %>% arrange(-n)
ccf_annotations %>% group_by(parcellation_division) %>% count() %>% arrange(-n) %>% as.data.frame()
ccf_annotations %>% group_by(parcellation_division, parcellation_category) %>% count() %>% arrange(-n) %>% as.data.frame()
ccf_annotations %>% group_by(parcellation_division, parcellation_division_color, parcellation_category) %>% count() %>% arrange(-n) %>% as.data.frame()

#group annotations with same color
ccf_annotations %<>% mutate(parcellation_division_remap = parcellation_division)
#VS is ventricular systems
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(parcellation_category == "VS", "VS", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(parcellation_category == "fiber tracts", "fiber tracts", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(parcellation_category == "brain-unassigned", "unassigned", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_color = if_else(parcellation_category == "brain-unassigned", "#000000", parcellation_division_color))

#manual mapping to match up with disscetion regions
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(grepl("^ACA", parcellation_structure), "ACA", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(grepl("^RSP", parcellation_structure), "RSP", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(grepl("^AUD", parcellation_structure), "AUD", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(grepl("^ENT", parcellation_structure), "ENT", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else("MOp" == parcellation_structure, "MOp", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else("SUB" == parcellation_structure, "SUB", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(grepl("^SSp", parcellation_structure), "S1", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(grepl("^VIS", parcellation_structure), "VIS", parcellation_division_remap))

ccf_annotations %<>% mutate(parcellation_division_remap = if_else("VISp" == parcellation_structure, "VISP", parcellation_division_remap))

ccf_annotations %<>% mutate(parcellation_division_remap = if_else("MY" == parcellation_division, "BS", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else("P" == parcellation_division, "BS", parcellation_division_remap))


ccf_annotations %<>% mutate(parcellation_division_remap = if_else(parcellation_division == "STR"& parcellation_structure %in% c("LSc","LSr","LSv","SF","SH"), 
                                                                  "LSX", parcellation_division_remap))
ccf_annotations %<>% mutate(parcellation_division_remap = if_else(parcellation_division == "STR"& parcellation_structure %in% c("CP"), 
                                                                  "STRd", parcellation_division_remap))


ccf_annotations %>% group_by(parcellation_division, parcellation_division_color, parcellation_category, parcellation_division_remap) %>% count() %>% arrange(-n) %>% as.data.frame()
ccf_annotations %>% group_by(parcellation_division_color, parcellation_division_remap) %>% count() %>% arrange(-n) %>% as.data.frame()
colnames(ccf_annotations)
joined %<>% left_join(ccf_annotations)
#write out for Figure_3_examine_MERFISH_recip_centroids.R - could be reduced in size by filtering NA's
joined %>% write_csv(paste0(base_folder, "spatial_per_cell_data.", top_hits %>% nrow(), "_pairs.csv.gz"))

joined %<>% mutate(matching_reciprocal = Zeng_merged_name == Macosko_merged_name)
joined %>% pull(matching_reciprocal) %>% sum(na.rm=T)

joined %>% filter(matching_reciprocal) %>% group_by(Macosko_merged_name) %>% count() %>% arrange(-n)

joined %>% group_by(is.na(parcellation_index), matching_reciprocal) %>% count()


#this is filtered on the Zeng side
slice_summary <- joined %>% group_by(slice) %>% summarize(total_cells = n(), 
                                                          Zeng_calls = sum(!is.na(Celltype_Zeng))/total_cells,
                                                          Mac_calls = sum(!is.na(Celltype_Macosko))/total_cells,
                                                          Zeng_recip_calls = sum(!is.na(Zeng_merged_name))/total_cells,
                                                          Mac_recip_calls = sum(!is.na(Macosko_merged_name))/total_cells,
                                                          matching_reciprocal = sum(matching_reciprocal, na.rm=T)/total_cells,
                                                          unique_merged_names = na.omit(union(Zeng_merged_name, Macosko_merged_name)) %>% unique() %>% length(), 
                                                          unique_clusters_all = na.omit(union(Celltype_Zeng, Celltype_Macosko)) %>% unique() %>% length(), 
                                                          proportion_merged_names = unique_merged_names/na.omit(union(joined$Zeng_merged_name, joined$Macosko_merged_name)) %>% unique() %>% length(),
                                                          proportion_clusters_seen = unique_clusters_all/na.omit(union(joined$Celltype_Zeng, joined$Celltype_Macosko)) %>% unique() %>% length()
)

cor.test(slice_summary$proportion_clusters_seen, slice_summary$proportion_merged_names) #0.96 correlation
cor.test(slice_summary$proportion_clusters_seen, slice_summary$matching_reciprocal, m='s') #0.5 correlation

#write out
slice_summary %>% select(slice, total_cells, proportion_Zeng_recip_calls = Zeng_recip_calls, proportion_Mac_recip_calls = Mac_recip_calls,  cells_with_matching_reciprocal_calls = matching_reciprocal, 
                         detected_reciprocal_clusters = unique_merged_names, proportion_of_all_reciprocals = proportion_merged_names)  %>% write_csv(paste0(base_folder, "spatial_slice_stats.", top_hits %>% nrow(), "_pairs.csv"))
correlation <- cor(slice_summary$unique_clusters_all, slice_summary$matching_reciprocal,m='s')

for_line_plot <- slice_summary %>% select(slice,matching_reciprocal,proportion_clusters_seen) %>% pivot_longer(cols=c(matching_reciprocal, proportion_clusters_seen))
for_line_plot %<>% mutate(name = if_else(name == "matching_reciprocal", "MERFISH cells mapped to reciprocal pairs","Clusters detected"))
for_line_plot %<>% mutate(slice_factor = factor(slice, levels = for_line_plot %>% pull(slice) %>% unique() %>% sort(decreasing = T)))
for_line_plot %<>% mutate(slice_numeric = as.numeric(slice_factor))
target_slices <- for_line_plot %>% filter(slice_factor %in% c(14, 24, 40, 60)) %>% pull(slice_factor) %>% unique()
target_slices_numeric <- for_line_plot %>% filter(slice_factor %in% c(14, 24, 40, 60)) %>% pull(slice_numeric) %>% unique()

pdf(paste0(base_folder, "spatial_slice_stats.", top_hits %>% nrow(), "_pairs.pdf"), height=5, width = 11)
ggplot(data= for_line_plot, aes(x=slice_numeric, y=value, color = name)) + 
  #geom_vline(xintercept = target_slices_numeric, color="darkgray") +
  geom_line() + geom_point() + 
  theme_bw() + theme(legend.position = "right") + labs(color="") + ylab("Proportion") + 
  xlab("Slice") + scale_x_continuous(labels = c("Anterior", "Posterior"), breaks = c(1,59)) +
  annotate("text",
           x = 1, 
           y = 0.6, 
           label = paste0("r: ", correlation %>% format(digits = 2)),
           hjust = 0, 
           vjust = 1
  )  + scale_color_manual(values=c("#795695", "#f2ba88")) 
dev.off()

#split by the two datasets
total_Zeng <- joined$Celltype_Zeng %>% unique() %>% length()
total_Mac <- joined$Celltype_Macosko %>% unique() %>% length()

#need to set some Zeng cell types to NA to match size
joined_cell_count_matched <- joined 
joined_cell_count_matched %<>% select(cell_label, slice, Celltype_Zeng, Celltype_Macosko)

IDs_to_remove_calls <- joined_cell_count_matched %>% filter(!is.na(Celltype_Zeng)) %>% sample_n(joined_cell_count_matched %>% filter(!is.na(Celltype_Macosko)) %>% nrow()) %>% pull(cell_label)
joined_cell_count_matched %<>% mutate(Celltype_Zeng = if_else(cell_label %in% IDs_to_remove_calls, Celltype_Zeng, NA))

joined_cell_count_matched %>% group_by(!is.na(Celltype_Zeng)) %>% count()
joined_cell_count_matched %>% group_by(!is.na(Celltype_Macosko)) %>% count()
#need to match the same number of cells for Zeng as Macosko
slice_summary_extra_counts <- joined_cell_count_matched %>% group_by(slice) %>% summarize(
  unique_clusters_Zeng = na.omit(Celltype_Zeng) %>% unique() %>% length(), 
  unique_clusters_Mac = na.omit(Celltype_Macosko) %>% unique() %>% length(),
  value = unique_clusters_Zeng/total_Zeng,
  unique_clusters_Mac_prop = unique_clusters_Mac/total_Mac, 
  detected_cluster_difference = unique_clusters_Zeng - unique_clusters_Mac
)
cor.test(slice_summary_extra_counts$unique_clusters_Mac, slice_summary_extra_counts$unique_clusters_Zeng,m='s')
#add in reciprocal best hits
slice_summary_extra_counts %<>% inner_join(slice_summary %>% select(slice, matching_reciprocal))
cor.test(slice_summary_extra_counts$detected_cluster_difference, slice_summary_extra_counts$matching_reciprocal, m='s')
cor.test(slice_summary_extra_counts$unique_clusters_Zeng, slice_summary_extra_counts$matching_reciprocal, m='s')
cor.test(slice_summary_extra_counts$unique_clusters_Mac, slice_summary_extra_counts$matching_reciprocal, m='s')
plot(slice_summary_extra_counts$detected_cluster_difference, slice_summary_extra_counts$matching_reciprocal, m='s')

correlation <- cor(slice_summary_extra_counts$detected_cluster_difference, slice_summary_extra_counts$matching_reciprocal, m='s')

pdf(paste0(base_folder, "spatial_cluster_count_vs_recips.pdf"), height=5, width = 5)
ggplot(data = slice_summary_extra_counts, aes(x = detected_cluster_difference, y = matching_reciprocal)) + 
  geom_point() + theme_bw() + 
  xlab("Diversity mismatch\n(excess clusters detected in cells vs. nuclei derived mappings)") + 
  ylab("Reproducable cell calls\n(proportion of cells with matching reciprocal status)") +
  annotate("text",
           x = 700, 
           y = 0.45, 
           label = paste0("r: ", correlation %>% format(digits = 2)),
           hjust = 0, 
           vjust = 1
  )
dev.off()


#it would be better if we had x y z for all cells - missing 394k that weren't mapped in Zeng
joined %>% filter(is.na(x)) %>% nrow()
for_plot <- joined %>% filter(!is.na(x)) %>% select(cell_label, x, y, Celltype_Macosko, Celltype_Zeng, Zeng_merged_name, Macosko_merged_name, slice, matching_reciprocal)
# Plot one reciprocal
for_plot %<>% filter(slice %in% c(18,19))
for_plot %<>% mutate(slice = paste0("Slice ", slice))

for_plot %>% filter(Zeng_merged_name == "Zeng|2670__Macosko|CholEx_Tcf24_Oxtr")

#find small set of cells with recip match
#CholEx_Tcf24_Oxtr
joined %>% filter(Celltype_Macosko == "Macosko|CholEx_Tcf24_Oxtr") %>% pull(Zeng_merged_name) %>% length()
joined %>% filter(Celltype_Macosko == "Macosko|CholEx_Tcf24_Oxtr") %>% group_by(Celltype_Zeng) %>% count()
joined %>% filter(Celltype_Macosko == "Macosko|CholEx_Tcf24_Oxtr") %>% group_by(slice) %>% count()
joined %>% filter(Zeng_merged_name == "Zeng|2670__Macosko|CholEx_Tcf24_Oxtr")%>% group_by(slice) %>% count()
joined %>% filter(Celltype_Zeng == "Zeng|2670") %>% group_by(Celltype_Zeng) %>% count()


pdf(paste0(base_folder, "spatial_CholEx_Tcf24_Oxtr.pdf"), height=5, width = 10)
ggplot(data = for_plot, aes(x=x, y=-1*y)) + geom_point(size=.03,  color ="gray90") + theme_classic() +    
  #geom_point(data = for_plot %>% filter(Celltype_Macosko == "Macosko|CholEx_Tcf24_Oxtr", Zeng_merged_name != "Zeng|2670__Macosko|CholEx_Tcf24_Oxtr"), size = 0.7, color="#E21E25") +  #c("#4450A2", "#E21E25"))
  geom_point(data = for_plot %>% filter(Zeng_merged_name == "Zeng|2670__Macosko|CholEx_Tcf24_Oxtr"), size = 0.9, color="black") + 
  geom_point(data = for_plot %>% filter(Celltype_Zeng == "Zeng|2670", is.na(Celltype_Macosko) | Celltype_Macosko != "Macosko|CholEx_Tcf24_Oxtr"), size = 0.9, color="#4450A2") +  #c("#4450A2", "#E21E25"))
  coord_equal() + theme(legend.position="none") +
  facet_wrap(. ~ slice) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(
    axis.line = element_line(color = 'white'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  ) 
dev.off()

#Figure 5 plotting three clusters that lined up with the BARseq data. The Macosko data is pulled manually from the website.
for_plot <- joined %>% filter(!is.na(x)) %>% select(cell_label, x, y, Celltype_Macosko, Celltype_Zeng, Zeng_merged_name, Macosko_merged_name, slice, matching_reciprocal)
target_clusters_for_pdf <- c("2866", "5231", "4218")
target_cluster <- "4218"
target_cluster <- paste0("Zeng|", target_cluster)
max_slice <- for_plot %>% filter(Celltype_Zeng == target_cluster) %>% group_by(slice) %>% count() %>% ungroup() %>% slice_max(n) %>% pull(slice)
if (target_cluster == "Zeng|4218") { max_slice <- 24 } #this slice best fits the barseq
if (target_cluster == "Zeng|5231") { max_slice <- 30 } #this slice best fits the barseq
pdf(paste0(base_folder, "spatial_", target_cluster,"_slice_", max_slice,".pdf"), height=5, width = 6)
for_plot <- joined %>% filter(!is.na(x)) %>% select(cell_label, x, y, Celltype_Macosko, Celltype_Zeng, Zeng_merged_name, Macosko_merged_name, slice, matching_reciprocal)
for_plot %<>% filter(slice == max_slice)
ggplot(data = for_plot, aes(x=x, y=-1*y)) + geom_point(size=.03,  color ="gray90") + theme_classic() +    
  geom_point(data = for_plot %>% filter(Celltype_Zeng == target_cluster), size = 0.15, color="#4450A2") + 
  coord_equal() + theme(legend.position="none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(
    axis.line = element_line(color = 'white'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  ) 
dev.off()

## plot the reciprocal hit cells
for_plot <- joined %>% filter(!is.na(x)) %>% select(cell_label, x, y, Zeng_merged_name, Macosko_merged_name, slice, matching_reciprocal, parcellation_division_remap, parcellation_division_color)
#convert na's to unassigned
for_plot %<>% mutate(parcellation_division_remap = if_else(is.na(parcellation_division_remap), "unassigned", parcellation_division_remap))
for_plot %<>% mutate(parcellation_division_color = if_else(is.na(parcellation_division_color), "#000000", parcellation_division_color))
set.seed(42)

for_plot %>% filter(slice == 40) %>% group_by(matching_reciprocal) %>% count()

#target_slices <- c(14, 24, 40, 60)
#target_slices <- c("46", '45', '43', '42'  ,"40", '39','38', "37")
for_plot %<>% filter(slice %in% target_slices)

for_plot %>% pull(slice) %>% unique() %>% length()
for_plot %<>% mutate(slice_factor = factor(slice, levels = for_plot %>% pull(slice) %>% unique() %>% sort(decreasing = T)))
for_plot %<>% mutate(matching_reciprocal = if_else(is.na(matching_reciprocal), FALSE, matching_reciprocal))


pdf(paste0(base_folder, "spatial_four_recips.", top_hits %>% nrow(), "_pairs.pdf"), height=7, width = 30)
ggplot(data = for_plot, aes(x = x, y = -1 * y, color = factor(matching_reciprocal))) +
  geom_point(size=.02) +
  geom_point(data = for_plot %>% filter(matching_reciprocal), size = 0.02, color="black") + 
  theme_classic() +
  coord_equal() +
  facet_wrap(. ~ slice_factor, nrow = 1) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(color = 'white'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  ) +
  scale_color_manual(
    values = c("gray", "black"),
    labels = c("Non-matching", "Matching"),
    name = "Reciprocal status",
    guide = guide_legend(override.aes = list(size = 5))
  )
dev.off()

#Coarse region atlas color
for_plot %<>% filter(!parcellation_division_remap %in% c("fiber tracts", "unassigned", "VS"))
region_colors <- for_plot %>% select(parcellation_division_remap, parcellation_division_color) %>% distinct() %>% as.data.frame()

#division colored, all points with data
pdf(paste0(base_folder, "MERFISH_colored_divsions.pdf"), height=7, width = 30)
ggplot(data = for_plot, aes(x = x, y = -1 * y, color = parcellation_division_color)) + 
  geom_point(size = 0.02) +
  scale_color_identity(
    name = "Structure",
    labels = region_colors$parcellation_division_remap,
    breaks = region_colors$parcellation_division_color,
    guide = guide_legend(override.aes = list(size = 5))
    #guide = 'legend'
  ) + 
  coord_equal() + 
  theme_classic() +
  facet_wrap(. ~ slice_factor, nrow = 1) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(color = 'white'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  ) 
dev.off()


################################
################################
#cut up into bins for heatmap view
################################
################################

for_plot <- joined %>% filter(!is.na(x)) %>% select(cell_label, x, y, slice, matching_reciprocal)
for_plot %>% pull(x) %>% range()
for_plot %>% pull(y) %>% range()
ratio <- (max(for_plot %>% pull(y)) - min(for_plot %>% pull(y)))/(max(for_plot %>% pull(x)) - min(for_plot %>% pull(x)))
x_bins = 110
y_bins = round(x_bins * ratio)

for_plot %<>% mutate(bin_x = cut_interval(x, n = x_bins))
for_plot %<>% mutate(bin_y = cut_interval(-1*y, n = y_bins))

for_plot %<>% group_by(bin_x, bin_y, slice) %>% summarize(n=n(), prop_recip = sum(matching_reciprocal, na.rm=T)/n())

for_plot %<>% filter(slice %in% target_slices)
for_plot %<>% mutate(slice_factor = factor(slice, levels = for_plot %>% pull(slice) %>% unique() %>% sort(decreasing = T)))

pdf(paste0(base_folder, "spatial_four_recips_proportions.", top_hits %>% nrow(), "_pairs.pdf"), height=7, width = 32)
ggplot(data=for_plot, aes(x=bin_x, y=bin_y, fill = prop_recip) ) +
  geom_tile() + scale_fill_continuous(type = "viridis", na.value='white') + theme_bw() +
  labs(fill = "P(reciprocal match)") + coord_equal(expand=FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(
    axis.line = element_line(color = 'white'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  ) + facet_wrap(.~slice_factor, nrow=1)
dev.off()

##########################
#test enrichment by coarse region
###########################
joined %<>% mutate(matching_reciprocal = if_else(is.na(matching_reciprocal), FALSE, matching_reciprocal))
whole_dataset <- joined %>% group_by(matching_reciprocal)  %>% count()
whole_dataset %<>% pivot_wider(names_from=matching_reciprocal, values_from=n, names_prefix = "dataset_recip_hit_")
by_region <- joined %>% mutate(parcellation_division_remap = if_else(is.na(parcellation_division_remap), "unassigned", parcellation_division_remap)) %>%
  group_by(parcellation_division_remap, matching_reciprocal)  %>% count()
by_region %<>% pivot_wider(names_from=matching_reciprocal, values_from=n, names_prefix = "recip_hit_")
by_region %<>% mutate(dataset_recip_hit_FALSE = whole_dataset %>% pull(dataset_recip_hit_FALSE))
by_region %<>% mutate(dataset_recip_hit_TRUE = whole_dataset %>% pull(dataset_recip_hit_TRUE))
by_region %<>% mutate(total_cells_in_region = recip_hit_TRUE + recip_hit_FALSE)
by_region %<>% mutate(percent_recip = recip_hit_TRUE / total_cells_in_region)
by_region %<>% ungroup() %>% mutate(dataset_total_cells = sum(total_cells_in_region))
by_region %<>% mutate(dataset_prop_recip = dataset_recip_hit_TRUE/ dataset_total_cells)
#duplicaed code from regional_enrichment_tests
by_region %<>% mutate(p_enriched = phyper(recip_hit_TRUE, dataset_recip_hit_TRUE, dataset_recip_hit_FALSE, total_cells_in_region, lower.tail = FALSE)  + 
                        dhyper(recip_hit_TRUE, dataset_recip_hit_TRUE, dataset_recip_hit_FALSE, total_cells_in_region))

by_region %<>% mutate(p_depleted = phyper(recip_hit_TRUE, dataset_recip_hit_TRUE, dataset_recip_hit_FALSE, total_cells_in_region, lower.tail = TRUE))
by_region %<>% mutate(difference_from_dataset_average = percent_recip-dataset_prop_recip)

by_region %>% as.data.frame()

by_region %<>% mutate(Dataset= "MERFISH")
by_region %<>% mutate(fold_change = percent_recip/dataset_prop_recip) #not used
by_region %<>% ungroup()
dir.create(here(base_folder, "top_hit_enrichments"))

#used for plotting in regional_enrichment_tests.R
by_region %>% write_csv(here(base_folder, "top_hit_enrichments", paste0("MERFISH_regional_enrichment_results.", top_hits %>% nrow(), "_pairs.csv")))


