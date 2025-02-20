library(MetaNeighbor)
library(MetaMarkers)
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
Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)

results_base_folder <- here("results")

markers_Macosko <- read_csv(file.path(results_base_folder, "MetaMarkers_Macosko_cluster/Macosko_markers.csv.gz"))

markers_Zeng <- read_csv(file.path(results_base_folder, "MetaMarkers_Zeng_cluster/Zeng_markers.csv.gz"))
markers_Zeng %<>% mutate(cell_type = as.character(cell_type))
markers_Zeng

#load in tophits for alignment
#cell_type needs to be mapped
base_folder <- here("results/full_run_ZengAWS.1698256033/") #updated Allen data - 2009 recips
tophits <- read_csv(file.path(base_folder, "top_hits.0.95.csv"))
tophits
tophits %<>% filter(Match_type == "Reciprocal_top_hit")
tophits %>% as.data.frame()
#check it's Zeng then Mac
tophits %>% filter(grepl("Mac", `Study_ID|Celltype_1`))
tophits %>% filter(grepl("Zeng", `Study_ID|Celltype_2`))
tophits %<>% mutate(Macosko_full = get_cell_type( `Study_ID|Celltype_2`))
tophits %<>% mutate(Zeng_full = get_cell_type( `Study_ID|Celltype_1`))
tophits %<>% mutate(joined_label = paste0(Zeng_full, "__", Macosko_full))

#filter for best recip matches
markers_Zeng %<>% inner_join(tophits %>% select(cell_type = Zeng_full, joined_label))
markers_Zeng %<>% select(-cell_type) %>% rename(cell_type = joined_label)

markers_Macosko %<>% inner_join(tophits %>% select(cell_type = Macosko_full, joined_label))
markers_Macosko %<>% select(-cell_type) %>% rename(cell_type = joined_label)

setdiff(markers_Zeng %>% pull(cell_type), markers_Macosko %>% pull(cell_type))
setdiff(markers_Macosko %>% pull(cell_type), markers_Zeng %>% pull(cell_type))
common_labels <- intersect(markers_Zeng %>% pull(cell_type), markers_Macosko %>% pull(cell_type))
length(common_labels) #number of reciprocal pairs

markers_Zeng %<>% filter(cell_type %in% common_labels)
markers_Macosko %<>% filter(cell_type %in% common_labels)

markers = list(
  zeng = markers_Zeng,
  macosko = markers_Macosko
)

meta_markers = make_meta_markers(markers, detailed_stats = TRUE)
#write out
dir.create(file.path(results_base_folder, paste0("MetaMarkers_", length(common_labels), "_reciprocals")))
write_csv(meta_markers, file = file.path(results_base_folder, paste0("MetaMarkers_", length(common_labels), "_reciprocals/MetaMarkers_default.csv.gz")))
meta_markers %>% filter(recurrence == 2) %>% group_by(cell_type) %>% count() %>% arrange(-n)
meta_markers %>% filter(recurrence == 2) %>% group_by(cell_type) %>% count() %>% arrange(n)
write_csv(meta_markers %>% filter(recurrence == 2), file = file.path(results_base_folder, paste0("MetaMarkers_", length(common_labels), "_reciprocals/MetaMarkers_default_recurrence_2.csv.gz")))

# pareto_markers = get_pareto_markers(meta_markers, "Zeng|1005__Macosko|Inh_Sox8_Cyp26b1", min_recurrence=2)
# plot_pareto_markers(meta_markers, "Zeng|1005__Macosko|Inh_Sox8_Cyp26b1", min_recurrence = 0)