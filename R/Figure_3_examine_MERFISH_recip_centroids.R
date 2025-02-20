library(MetaNeighbor) #for heatmap code
library(stringr)
library(reshape2)
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
Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)

base_folder <- here("results/full_run_ZengAWS.1698256033/") #updated Allen data - 2009 recips

all_joined_centroids <- tibble()
all_joined_calls <- tibble()
#set to also run on the more strict 612 top hits
for (file in list.files(base_folder, pattern = "spatial_per_cell_data")) {
  print(file)
  joined <- read_csv(here(base_folder, file))
  all_joined_calls %<>% bind_rows(joined %>% mutate(source = file))
  
  centroids_Mac <- joined %>% group_by(Macosko_merged_name) %>% summarize(x = mean(x, na.rm=T), y = mean(y, na.rm=T), z = mean(z, na.rm=T))
  centroids_Zeng <- joined %>% group_by(Zeng_merged_name) %>% summarize(x = mean(x, na.rm=T), y = mean(y, na.rm=T), z = mean(z, na.rm=T))
  
  joined_centroids <- inner_join(centroids_Mac, centroids_Zeng, by = join_by(Macosko_merged_name == Zeng_merged_name), suffix = c(".Mac", ".Zeng"))
  
  #could use some refactoring since I also create the full all v all distance matrix
  joined_centroids %<>% rowwise() %>%
    mutate(distance = sqrt((x.Mac - x.Zeng)^2 + (y.Mac - y.Zeng)^2 + (z.Mac - z.Zeng)^2)) 
  joined_centroids %<>% mutate(source = file)
  all_joined_centroids %<>% bind_rows(joined_centroids)
}
print(file)
#some data duplication above, use the last file loaded 
for_all_centroids <- all_joined_calls %>% filter(source == file) %>% filter(!is.na(x_ccf))
all_centroids <- bind_rows(for_all_centroids %>% group_by(cluster_id = Celltype_Zeng) %>% summarize(x = mean(x, na.rm=T), y = mean(y, na.rm=T), z = mean(z, na.rm=T)),
                           for_all_centroids %>% group_by(cluster_id = Celltype_Macosko) %>% summarize(x = mean(x, na.rm=T), y = mean(y, na.rm=T), z = mean(z, na.rm=T)))
all_centroids %>% filter(!is.na(cluster_id))

filtered_centroids <- all_centroids %>%
  select(x, y, z)

# Compute the Euclidean distance matrix
distance_matrix <- dist(all_centroids %>% select(x, y, z), method = "euclidean")

distance_matrix <- as.matrix(distance_matrix)
dim(distance_matrix)

colnames(distance_matrix) <- all_centroids %>% pull(cluster_id)
rownames(distance_matrix) <- all_centroids %>% pull(cluster_id)
distance_matrix_raw <- distance_matrix
#scale between 0-1
distance_matrix <- distance_matrix/max(distance_matrix)
#invert so it's similarity
distance_matrix <- 1-distance_matrix

hist(distance_matrix) #normalized to zero to 1

all_joined_centroids %>% pull(distance) %>% hist()
all_joined_centroids %>% pull(distance) %>% range(na.rm=T)
#for statistics in manuscript
all_joined_centroids %>% group_by(source) %>% summarize(mean=mean(distance,na.rm=T), median = median(distance,na.rm=T), n=n(), na = sum(is.na(distance)))


summary_stats <- all_joined_calls %>% group_by(source) %>% summarize(
  Mac_calls = sum(!is.na(Celltype_Macosko)),
  Mac_recip_calls = sum(!is.na(Macosko_merged_name)),
  Zeng_calls = sum(!is.na(Celltype_Zeng)),
  Zeng_recip_calls = sum(!is.na(Zeng_merged_name)),
  Mac_proportion = Mac_recip_calls/Mac_calls,
  Zeng_proportion = Zeng_recip_calls/Zeng_calls
)
summary_stats %<>% mutate(n_clusters = as.numeric(str_extract(source, "\\d+")))
summary_stats %<>% mutate(Zeng_odds = Zeng_proportion/(n_clusters/5322))
summary_stats %<>% mutate(Macosko_odds = Mac_proportion/(n_clusters/5030))
summary_stats

#perform cluster size correlations in recips and non recips
Mac_call_counts <- all_joined_calls %>% group_by(source, Celltype_Macosko) %>% count()
Zeng_call_counts <- all_joined_calls %>% group_by(source, Celltype_Zeng) %>% count()

#add in size info to compare
merged_obs <- read_csv(paste0( base_folder, "merged.obs.csv.zip"))
merged_obs %<>% rename(cell_id = "...1")
cell_type_sizes <- merged_obs %>% group_by(cell.type, study_id) %>% count()
cell_type_sizes %<>% ungroup()
cell_type_sizes %<>% mutate(cell.type = paste0(study_id, "|", cell.type))

Mac_call_counts <- inner_join(Mac_call_counts, cell_type_sizes %>% mutate(Celltype_Macosko = cell.type), by = "Celltype_Macosko", suffix = c(".MERFISH", ".original"))
Mac_call_counts %>% group_by(source) %>% summarize(cor = cor(n.MERFISH, n.original, m='s'))
#split by recip status
Mac_call_counts %<>% left_join(all_joined_calls %>% select(source, Celltype_Macosko, Macosko_merged_name) %>% distinct() %>% filter(!is.na(Macosko_merged_name)))
Mac_call_counts %>% group_by(source, is_recip=!is.na(Macosko_merged_name)) %>% summarize(cor = cor(n.MERFISH, n.original, m='s'), median_n.original = median(n.original), median_n.MERFISH = median(n.MERFISH))

Zeng_call_counts <- inner_join(Zeng_call_counts, cell_type_sizes %>% mutate(Celltype_Zeng = cell.type), by = "Celltype_Zeng", suffix = c(".MERFISH", ".original"))
Zeng_call_counts %>% group_by(source) %>% summarize(cor = cor(n.MERFISH, n.original, m='s'))
Zeng_call_counts %<>% left_join(all_joined_calls %>% select(source, Celltype_Zeng, Zeng_merged_name) %>% distinct() %>% filter(!is.na(Zeng_merged_name)))
Zeng_call_counts %>% group_by(source, is_recip=!is.na(Zeng_merged_name)) %>% summarize(cor = cor(n.MERFISH, n.original, m='s'), median_n.original = median(n.original), median_n.MERFISH = median(n.MERFISH))

all_joined_centroids %<>% mutate(source = gsub("spatial_per_cell_data.", "", source))
all_joined_centroids %<>% mutate(source = gsub(".csv.gz", "", source))
all_joined_centroids %>% group_by(source) %>%  summarize(n=n(), median_centroid_distance = median(distance, na.rm = T)) %>% arrange(-median_centroid_distance)
all_joined_centroids %<>% mutate(source = factor(source, levels = c("972_pairs", "2009_pairs", "612_pairs")))
#remove for simplicity
all_joined_centroids %<>% filter(source != "612_pairs")

pdf(paste0(base_folder, "spatial_centroid_distances.pdf"), width=6, height = 4)
ggplot(data = all_joined_centroids, aes(fill = source, x = distance)) + 
  geom_density() + theme_bw() + xlab("Distance between centroids of the pairs") + facet_wrap(. ~ source, nrow=2) +
  geom_vline(data = means, aes(xintercept = mean_distance), linetype = "dashed", color = "black", size = 0.8) +
  scale_x_continuous(limits = c(0, NA), expand = c(0, 0))
dev.off()



one_vs_one_AUROCs <- read_csv(paste0(base_folder, "aurocs_1v1.csv.gz"), guess_max = 200000000)
one_vs_one_AUROCs %<>% rename(target = "...1")
#melting to pull out more details
one_vs_one_AUROCs_melted <- melt(one_vs_one_AUROCs, na.rm = T, value.name= "auroc", variable.name = "reference") %>% tibble()
one_vs_one_AUROCs_melted %<>% mutate(reference = as.character(reference), target = as.character(target))
one_vs_one_AUROCs_melted %<>% group_by(reference) %>% slice_min(auroc) %>% select(reference, second_best = target)
one_vs_one_AUROCs_melted %<>% filter(reference %in% rownames(distance_matrix))
one_vs_one_AUROCs_melted %<>% filter(second_best %in% rownames(distance_matrix))

one_vs_one_AUROCs_melted %<>% rowwise() %>% mutate(distance = distance_matrix_raw[reference, second_best]) #don't use the normalized distance matrix

all_joined_centroids_for_1v1 <- all_joined_centroids
all_joined_centroids_for_1v1 %<>% separate(Macosko_merged_name, into = c("Zeng_cluster", "Macosko_cluster"), sep = "__")
all_joined_centroids_for_1v1 %<>% inner_join(one_vs_one_AUROCs_melted %>% select(Zeng_cluster = reference, distance_Zeng_reference = distance))
all_joined_centroids_for_1v1 %<>% inner_join(one_vs_one_AUROCs_melted %>% select(Macosko_cluster = reference, distance_Macosko_reference = distance))
all_joined_centroids_for_1v1 %>% group_by(distance < distance_Macosko_reference) %>% count()
all_joined_centroids_for_1v1 %>% group_by(distance < distance_Zeng_reference) %>% count()

for_plot <- bind_rows(all_joined_centroids_for_1v1 %>% select(Zeng_cluster, Macosko_cluster, source, distance) %>% mutate(type = "Original"),
                      all_joined_centroids_for_1v1 %>% select(Zeng_cluster, Macosko_cluster, source, distance = distance_Zeng_reference) %>% mutate(type = "Next, Zeng reference"),
                      all_joined_centroids_for_1v1 %>% select(Zeng_cluster, Macosko_cluster, source, distance = distance_Macosko_reference) %>% mutate(type = "Next, Macosko reference"))

for_plot %>% filter(type == "Next, Macosko reference", source == "2009_pairs") %>% pull(distance) %>% median()
for_plot %>% filter(type == "Next, Zeng reference", source == "2009_pairs") %>% pull(distance) %>% median()

for_plot %>% filter(type == "Next, Macosko reference", source == "2009_pairs") %>% pull(distance) %>% mean()
for_plot %>% filter(type == "Next, Zeng reference", source == "2009_pairs") %>% pull(distance) %>% mean()

wilcox.test(for_plot %>% filter(type == "Original", source == "2009_pairs") %>% pull(distance), 
            for_plot %>% filter(type == "Next, Macosko reference", source == "2009_pairs") %>% pull(distance))
t.test(for_plot %>% filter(type == "Original", source == "2009_pairs") %>% pull(distance), 
       for_plot %>% filter(type == "Next, Macosko reference", source == "2009_pairs") %>% pull(distance))
