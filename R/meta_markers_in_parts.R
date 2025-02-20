library(here)
library(tidyr)
library(stringr)
library(org.Mm.eg.db)
library(MetaMarkers)
library(magrittr)
library(dplyr)
library(zellkonverter)
library(SingleCellExperiment)
library(readr)
library(conflicted)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)

#this loads up the h5ad files split by gene to reduce the memory load for the metamarkers analysis, still takes up a lot of memory
results_base_folder <- here("results")

###########################
#Set target here. 
#target_dataset = "Macosko"
target_dataset = "Zeng"
###########################


if (target_dataset == "Macosko") {
  name_of_run <- "MetaMarkers_Macosko_cluster"
  subset_path <- here("data/whole_mouse_brain/processed/macosko/subsets")
  target_pattern = "Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.cpm.part"
} else if (target_dataset == "Zeng") {
  name_of_run <- "MetaMarkers_Zeng_cluster"
  subset_path <- here("data/whole_mouse_brain/processed/zeng/subsets")
  target_pattern = "AIT21.0.merged.with_multiome.cpm.part"
}

dir.create(file.path(results_base_folder, name_of_run))

file_list <- list.files(subset_path, pattern = target_pattern, full.names = TRUE)
file_list <- sort(file_list)

markers_all <- tibble()

#this will accumulate memory in spite of the gc calls, so it may need to be ran across multiple sessions
for(part_file in file_list) {
  print(part_file)
  part_num <- str_extract(part_file, "part_\\d+")
  print(part_num)
  h5ad_part <- NULL
  gc()
  h5ad_part <- readH5AD(part_file, verbose=TRUE)
  print("Done reading")
  
  #for testing
  #h5ad_part <- h5ad_part[, sample(ncol(h5ad_part), round(ncol(h5ad_part) * 0.01))]
  
  if (target_dataset == "Macosko") {
    markers_part <- compute_markers(assay(h5ad_part, "X"), h5ad_part$ClusterNm)
  } else if (target_dataset == "Zeng") {
    markers_part <- compute_markers(assay(h5ad_part, "X"), h5ad_part$cl)  
  }
  
  print("Done compute markers")
  markers_part %<>% mutate(filename = part_file)
  
  #write out part in case it is interrupted due to memory
  write_csv(markers_part, path = file.path(results_base_folder, name_of_run, paste0(target_dataset, "_markers.", part_num, ".csv.gz")))
  
  markers_all <- bind_rows(markers_all, markers_part)
  rm(h5ad_part, markers_part)
  gc()
}

write_csv(markers_all, path = file.path(results_base_folder, name_of_run, paste0(target_dataset, "_markers.csv.gz")))

#code for gene_symbols
markers_all %<>% select(-group)
mapping <- mapIds(org.Mm.eg.db, markers_all %>% pull(gene) %>% unique(), 'SYMBOL', 'ENSEMBL')

mapping <- tibble(gene = names(mapping), gene_symbol = mapping) %>% distinct()
mapping %>% group_by(gene) %>% count() %>% filter(n != 1)
mapping %>% group_by(gene_symbol) %>% count() %>% filter(n != 1)

markers_all %<>% left_join(mapping)

markers_all %<>% select(-filename)
markers_all %<>% select(gene_symbol, everything())

write_csv(markers_all, path = file.path(results_base_folder, name_of_run, paste0(target_dataset, "_markers_with_gene_symbols.csv.gz")))
