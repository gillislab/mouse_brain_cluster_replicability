library(MetaNeighbor)
library(SingleCellExperiment)
library(magrittr)
library(dplyr)
library(S4Vectors)
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


base_data_folder <- here()

#This code is for the 109 gene pretrained models, the larger versions based on highly variable genes are at:
#Zeng_pretrained <- here('data/whole_mouse_brain/processed/zeng/merged_Zeng_AWS.Oct2023.pretrained_Zeng.csv.gz')
#Macosko_pretrained <- here('data/whole_mouse_brain/processed/macosko/merged_Zeng_AWS.Oct2023.pretrained_Macosko.csv.gz')

base_barseq_folder <- paste0(base_data_folder, "/data/whole_mouse_brain/processed/barseq/")

result_folder <- here("results", paste0("barseq_109_pretrained_R.", round(as.numeric(Sys.time()))))
dir.create(result_folder)

#load labels from https://data.mendeley.com/datasets/8bhhk7c5n9/1
labels <- read_csv(paste0(base_barseq_folder, "labels_20211208_withgaba_nonexc.csv"))
labels %>% filter(cluster == "Sst")
labels %>% pull(cluster) %>% unique() %>% sort() %>% length() #127 after removing NA value

#propigate labels down to cluster level for NA's
labels %<>% mutate(subclass = if_else(is.na(subclass), class, subclass))
labels %<>% mutate(cluster = if_else(is.na(cluster), subclass, cluster))
labels %>% pull(cluster) %>% unique() %>% sort()
labels %>% pull(cluster) %>% unique() %>% sort() %>% length()
labels %>% filter(cluster == "Sst")

#load rds from /data/fischer/ct_assessment/
bar_seq_data <- readRDS(paste0(base_barseq_folder, "barseq_210630.rds"))

rownames(bar_seq_data)
#remove unnamed genes
bar_seq_data <- bar_seq_data[rownames(bar_seq_data) != "unused-1", ]
bar_seq_data <- bar_seq_data[rownames(bar_seq_data) != "unused-2", ]
bar_seq_data <- bar_seq_data[rownames(bar_seq_data) != "unused-3", ]
bar_seq_data <- bar_seq_data[rownames(bar_seq_data) != "unused-4", ]
bar_seq_data <- bar_seq_data[rownames(bar_seq_data) != "unused-5", ]

#join labels
col_data <- colData(bar_seq_data) %>% as.data.frame() %>% rownames_to_column(var = "row_name") %>% tibble()
col_data %<>% mutate(sample = paste0(slice, "_", row_name))
col_data %<>% select(sample, everything())
col_data %<>% left_join(labels %>% select(-notes))
col_data %>% group_by(is.na(cluster)) %>% count()
col_data %>% group_by(is.na(cluster), is.na(ccf_x)) %>% count() #We don't filter on NA CCF coords (Cross expression does)

top_region_barseq <- col_data %>% filter(!is.na(cluster), !is.na(ccf_x)) %>% group_by(cluster, barseq_region_top = ccf_name) %>% 
  count() %>% group_by(cluster) %>% slice_max( n, with_ties = FALSE)
top_region_barseq %>% write_csv(paste0(base_barseq_folder, "top_region_per_cluster.csv"))


# Convert tibble to DataFrame
col_data_df <- DataFrame(col_data)
rownames(col_data_df) <- col_data_df$row_name
colData(bar_seq_data) <- col_data_df

bar_seq_data <- bar_seq_data[, !is.na(colData(bar_seq_data)$cluster)]
dim(bar_seq_data)

#check counts, slightly different due to the unused gene removal
counts <- tibble(gene_count = (counts(bar_seq_data)>0) %>% colSums(), read_count = counts(bar_seq_data) %>% colSums())
min(counts$gene_count)
min(counts$read_count)

#convert gene symbols
gene_map <- read_csv(paste0(base_barseq_folder, "barseq_gene_ids.csv"))
setdiff(rownames(bar_seq_data), gene_map$names)
setdiff(gene_map$names, rownames(bar_seq_data))
ensmbl_ordered <- tibble(names = rownames(bar_seq_data)) %>% left_join(gene_map) %>% pull(ENSEMBL)
rownames(bar_seq_data) <- ensmbl_ordered
dim(bar_seq_data)

#load pretrained model generated in python
ptrained_Zeng_109 = read_csv(paste0(base_data_folder, '/data/whole_mouse_brain/processed/zeng/subsets/AIT21.0.merged.with_multiome.pretrained_Zeng_109_genes.csv.gz'))
#check symbols
intersect(ptrained_Zeng_109 %>% pull('...1'), ensmbl_ordered) %>% length()
ptrained_Zeng_109 %<>%
  as.data.frame() %>%
  column_to_rownames('...1') 

ptrained_Mac_109 = read_csv(paste0(base_data_folder, "/data/whole_mouse_brain/processed/macosko/subsets/Macosko_Mouse_Atlas_Single_Nuclei.pretrained_Macosko_109_genes.csv.gz"))
#check gene count 
intersect(ptrained_Mac_109 %>% pull('...1'), ensmbl_ordered) %>% length()
ptrained_Mac_109 %<>%
  as.data.frame() %>%
  column_to_rownames('...1') 

bar_seq_data$study_id <- "BARseq"
#run pretrain
global_aurocs_Zeng <- MetaNeighborUS(
  trained_model = ptrained_Zeng_109, dat = bar_seq_data,
  study_id = bar_seq_data$study_id, cell_type = bar_seq_data$cluster,
  fast_version = TRUE
)

global_aurocs_Macosko <- MetaNeighborUS(
  trained_model = ptrained_Mac_109, dat = bar_seq_data,
  study_id = bar_seq_data$study_id, cell_type = bar_seq_data$cluster,
  fast_version = TRUE
)

dim(global_aurocs_Zeng)
dim(global_aurocs_Macosko)

write.table(global_aurocs_Zeng, gzfile(paste0(result_folder, "/Zeng_pretrained_aurocs.csv.gz")), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

write.table(global_aurocs_Macosko, gzfile(paste0(result_folder, "/Macosko_pretrained_aurocs.csv.gz")), sep = ",", row.names = TRUE, col.names = NA, quote = FALSE)

#write out cell labels used, for counting sizes
colData(bar_seq_data) %>% as.data.frame() %>% tibble() %>% group_by(cluster) %>% count() %>% write_csv(paste0(result_folder, "/BARseq_cluster_counts.csv"))
result_folder
