library(ggplot2)
library(data.table)
library(matrixStats)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)
library(readr)
library(conflicted)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)

#warning: duplicated loading code with coordinated expression analyses
subset_path_cell <- here("data", "whole_mouse_brain",  "processed", "zeng", "subsets/")
subset_path_nuc <- here("data", "whole_mouse_brain", "processed", "macosko", "subsets/")

#2009 recip hits
path_to_recip_hits <- here("results","full_run_ZengAWS.1698256033/")
top_hits_filename <- "top_hits.0.95.csv"

remove_constant_columns <- function(matrix_data) {
  # Identify constant columns
  constant_cols <- apply(matrix_data, 2, function(col) all(col == col[1]))
  
  matrix_without_constants <- matrix_data[, !constant_cols, drop = FALSE]
  
  return(matrix_without_constants)
}

#Pearson correlation between two matrices. Mirroring the R cor function, it performs column versus column correlation 
fast_pearson_matrix_cor <- function(matrixA, matrixB) {
  matrixA <- sweep(matrixA, 2, colMeans(matrixA), FUN = "-") 
  matrixB <- sweep(matrixB, 2, colMeans(matrixB), FUN = "-") 
  
  matrixA <- sweep(matrixA, 2, sqrt(colSums(matrixA^2)), FUN = "/") 
  matrixB <- sweep(matrixB, 2, sqrt(colSums(matrixB^2)), FUN = "/") 
  
  cr = crossprod(matrixA, matrixB)
  #dimnames(cr) = list(colnames(matrixA), colnames(matrixB))
  cr
}
top_hits <- read_csv(paste0(path_to_recip_hits, top_hits_filename))
top_hits %<>% select(-`...1`)
top_hits %<>% rowwise() %>% mutate(first = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[2],
                                   second = sort(c(`Study_ID|Celltype_1`, `Study_ID|Celltype_2`))[1])
top_hits %<>% mutate(Zeng_ID = first, Macosko_ID = second)
top_hits %<>% select(Zeng_ID, Macosko_ID, Match_type)
top_hits %>% group_by(Match_type) %>% count()
top_hits %<>% filter(Match_type == "Reciprocal_top_hit")
top_hits %<>% mutate(joined_ID = paste0(Zeng_ID, "__", Macosko_ID))
top_hits %<>% mutate(Zeng_ID = gsub(".*[|]","", Zeng_ID))
top_hits %<>% mutate(Macosko_ID = gsub(".*[|]","", Macosko_ID))

#some code duplication with other centroid code
single_cell_centroids <- read_csv(paste0(subset_path_cell, "AIT21.0.merged.with_multiome.log2_cpm.mean_centroids.csv.gz"))
single_cell_centroids[1:5,1:5]
dim(single_cell_centroids)

#############
#convert to a matrix
single_cell_centroids_matrix <- as.matrix(single_cell_centroids %>% select(-gene_id))
rownames(single_cell_centroids_matrix) <- single_cell_centroids %>% pull(gene_id)
dim(single_cell_centroids_matrix)
single_cell_centroids_matrix_full <- single_cell_centroids_matrix
#map names based on the recent run
single_cell_cols <- tibble(Zeng_ID = colnames(single_cell_centroids_matrix))
single_cell_cols %<>% mutate(row_number = row_number())
single_cell_cols %<>% inner_join(top_hits %>% select(Zeng_ID, joined_ID))
single_cell_cols %<>% arrange(row_number) #seems unnessicary but good to be certain

single_cell_centroids_matrix <- single_cell_centroids_matrix[, single_cell_cols %>% pull(Zeng_ID)]
colnames(single_cell_centroids_matrix) <- single_cell_cols %>% pull(joined_ID)
single_cell_centroids_matrix[1:5,1:5]
dim(single_cell_centroids_matrix)

######################
#do the same for single nuc
single_nuc_centroids <- read_csv(paste0(subset_path_nuc, "Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.log2_cpm.mean_centroids.csv.gz"))
single_nuc_centroids_matrix <- as.matrix(single_nuc_centroids %>% select(-gene_id))
rownames(single_nuc_centroids_matrix) <- single_nuc_centroids %>% pull(gene_id)
intersect(rownames(single_nuc_centroids_matrix), single_cell_centroids %>% pull(gene_id)) %>% length() #21649
single_nuc_centroids_matrix_full <- single_nuc_centroids_matrix #full version before 2009 set

#map names based on the recent run
single_nuc_cols <- tibble(Macosko_ID = colnames(single_nuc_centroids_matrix))
single_nuc_cols %<>% mutate(row_number = row_number())
single_nuc_cols %<>% inner_join(top_hits %>% select(Macosko_ID, joined_ID))
single_nuc_cols %<>% arrange(row_number) #seems unnessicary but good to be certain

single_nuc_centroids_matrix <- single_nuc_centroids_matrix[, single_nuc_cols %>% pull(Macosko_ID)]
colnames(single_nuc_centroids_matrix) <- single_nuc_cols %>% pull(joined_ID)

# Filter the matrix based on the non-zero rows
non_zero_rows <- apply(single_cell_centroids_matrix, 1, function(row) any(row != 0))
single_cell_centroids_matrix <- single_cell_centroids_matrix[non_zero_rows, ]

# Filter the matrix based on the non-zero rows
non_zero_rows <- apply(single_nuc_centroids_matrix, 1, function(row) any(row != 0))
single_nuc_centroids_matrix <- single_nuc_centroids_matrix[non_zero_rows, ]
dim(single_cell_centroids_matrix )
dim(single_nuc_centroids_matrix)

#set same order 
common_genes <- intersect(rownames(single_nuc_centroids_matrix), rownames(single_cell_centroids_matrix))
common_cell_types <- intersect(colnames(single_nuc_centroids_matrix), colnames(single_cell_centroids_matrix))
length(common_cell_types)
length(common_genes)

#common_cell_types <- head(common_cell_types) #uncomment for testing

single_nuc_centroids_matrix_full <- single_nuc_centroids_matrix_full[common_genes, ]
single_cell_centroids_matrix_full <- single_cell_centroids_matrix_full[common_genes, ]

single_nuc_centroids_matrix <- single_nuc_centroids_matrix[common_genes, common_cell_types]
single_cell_centroids_matrix <- single_cell_centroids_matrix[common_genes, common_cell_types]
dim(single_cell_centroids_matrix )
dim(single_nuc_centroids_matrix)

dim(single_nuc_centroids_matrix_full)
dim(single_cell_centroids_matrix_full)


################################
#transpose so we correlate genes, if statement so we don't do this multiple times
dim(single_nuc_centroids_matrix)

if (ncol(single_nuc_centroids_matrix) == nrow(top_hits)) {
  single_nuc_centroids_matrix <- t(single_nuc_centroids_matrix)
  single_cell_centroids_matrix <- t(single_cell_centroids_matrix)
}

#sparsity information
sum(single_cell_centroids_matrix == 0)/(length(single_cell_centroids_matrix))
sum(single_nuc_centroids_matrix == 0)/(length(single_nuc_centroids_matrix))

#get HVG info
hvgs <- read_csv(paste0(path_to_recip_hits, "merged.var.csv.zip"))
hvgs <- hvgs %>% pull("...1")
not_hvgs <- setdiff(colnames(single_cell_centroids_matrix), hvgs)
length(hvgs)
length(not_hvgs)

#########################
#filter out HVG's
single_cell_centroids_matrix_ranked <- apply(single_cell_centroids_matrix[,not_hvgs], 2, frank, ties.method = "average")
single_nuc_centroids_matrix_ranked <- apply(single_nuc_centroids_matrix[,not_hvgs], 2, frank, ties.method = "average")
dim(single_nuc_centroids_matrix_ranked)
dim(single_cell_centroids_matrix_ranked)

start_time <- Sys.time()
cor_matrix <- fast_pearson_matrix_cor(single_cell_centroids_matrix_ranked, single_nuc_centroids_matrix_ranked)
print(paste("Correlation elapsed time:", Sys.time() - start_time))

mean(diag(cor_matrix))
mean(cor_matrix)

cor_matrix[1:5,1:5] #not symmetric - need to do the ranking both directions or average the correlations
average_cor_matrix <- (cor_matrix + t(cor_matrix))/2
average_cor_matrix[1:5,1:5]
cor_matrix_ranked <- apply(average_cor_matrix, 2, frank, ties.method = "average")
mean(diag(cor_matrix_ranked))
#percentile of same-gene
mean(diag(cor_matrix_ranked))/nrow(cor_matrix_ranked)
mean(cor_matrix_ranked)
max(cor_matrix_ranked)

ranked_correlations <- tibble(gene_id = rownames(cor_matrix),  cor = diag(cor_matrix_ranked))
ranked_correlations %<>% arrange(-cor)
ranked_correlations %>% tail(20)
ranked_correlations %>% filter(cor == max(cor)) %>% nrow()
#their coordinated expression was higher than all correlations between that gene and any other
ranked_correlations %>% filter(cor == max(cor)) %>% nrow()/ranked_correlations %>% nrow()

ranked_correlations_to_write_out <- ranked_correlations
ranked_correlations_to_write_out %<>% rename(ranked_cor = cor)
ranked_correlations_to_write_out %<>% mutate(is_HVG = gene_id %in% hvgs)
ranked_correlations_to_write_out %<>% mutate(is_same_gene_best = ranked_cor == ranked_correlations %>% nrow())
ranked_correlations_to_write_out %>% write_csv(paste0(path_to_recip_hits, "coordinated_expression_by_gene.csv"))

###############
###############
#cut down the number of cell types and get
# mean spearman
# genes with top-1 hit
# average precentile
# exclude the HVG's?
all_results <- tibble()
#aborted_results <- all_results
set.seed(1)
target_sizes <- seq(nrow(single_cell_centroids_matrix), 1, by=-100)
target_sizes <- target_sizes[target_sizes != 709]
target_sizes <- target_sizes[target_sizes != 909]
target_sizes <- target_sizes[target_sizes != 1109]
target_sizes <- target_sizes[target_sizes != 1209]
target_sizes <- target_sizes[target_sizes != 1309]
target_sizes <- target_sizes[target_sizes != 1409]
target_sizes <- target_sizes[target_sizes != 1609]
target_sizes <- target_sizes[target_sizes != 1709]
target_sizes <- target_sizes[target_sizes != 1809]
target_sizes <- target_sizes[target_sizes != 1909]
target_sizes

for(cluster_count in target_sizes) {
  print(cluster_count)
  for(iteration in 1:20) {
    print(iteration)
    clusters_to_use <- sort(sample(rownames(single_cell_centroids_matrix), cluster_count))
    
    single_nuc_centroids_matrix_sampled <- single_nuc_centroids_matrix[clusters_to_use, ]
    single_cell_centroids_matrix_sampled <- single_cell_centroids_matrix[clusters_to_use, ]
    
    # Filter the matrix for constant genes
    single_nuc_centroids_matrix_sampled %<>% remove_constant_columns()
    single_cell_centroids_matrix_sampled %<>% remove_constant_columns()
    
    #set same order 
    common_genes <- intersect(colnames(single_cell_centroids_matrix_sampled), colnames(single_nuc_centroids_matrix_sampled))
    
    #### REMOVE HVG's here
    print("Removing HVGs")
    common_genes <- setdiff(common_genes, hvgs)
    
    single_cell_centroids_matrix_sampled <- single_cell_centroids_matrix_sampled[,common_genes]
    single_nuc_centroids_matrix_sampled <- single_nuc_centroids_matrix_sampled[,common_genes]
    
    single_cell_centroids_matrix_ranked <- apply(single_cell_centroids_matrix_sampled, 2, frank, ties.method = "average")
    single_nuc_centroids_matrix_ranked <- apply(single_nuc_centroids_matrix_sampled, 2, frank, ties.method = "average")
    
    cor_matrix <- fast_pearson_matrix_cor(single_cell_centroids_matrix_ranked, single_nuc_centroids_matrix_ranked)

    mean_diag_cor <- mean(diag(cor_matrix), na.rm=T)
    mean_whole_matrix <- mean(cor_matrix, na.rm=T)
    
    #symmetric version
    average_cor_matrix <- (cor_matrix + t(cor_matrix))/2
    
    cor_matrix_ranked <- apply(average_cor_matrix, 2, frank, ties.method = "average")
    #percentile of same-gene
    mean_percentile_same_gene <- mean(diag(cor_matrix_ranked), na.rm=T)/nrow(cor_matrix_ranked)
    
    ranked_correlations <- tibble(gene_id = rownames(cor_matrix_ranked),  cor = diag(cor_matrix_ranked), cor_abs = diag(cor_matrix))
    ranked_correlations %<>% filter(!is.na(cor_abs))
    
    top_1_proportion <- (ranked_correlations %>% filter(cor == max(cor)) %>% nrow())/(ranked_correlations %>% nrow())
    
    single_row <- tibble(iteration, cluster_count, mean_diag_cor, mean_whole_matrix, mean_percentile_same_gene, top_1_proportion, gene_count = ncol(single_nuc_centroids_matrix_ranked))
    all_results %<>% bind_rows(single_row)
    
    if (cluster_count ==  nrow(single_cell_centroids_matrix)) {
      print("Just doing one iteration of the full matrix")
      break
    }
  }
}

tail(all_results)
head(all_results)

#filter it
target_sizes


all_results_long <- all_results %>% pivot_longer(cols = mean_diag_cor:gene_count)


ggplot(data = all_results_long, aes(y=value, x=cluster_count)) + geom_point() + scale_x_reverse() + theme_bw() + facet_wrap(. ~ name, scales = "free_y")

all_results_long %<>% mutate(cluster_count_factor = factor(cluster_count, levels = all_results_long %>% pull(cluster_count) %>% unique()))
ggplot(data = all_results_long, aes(y=value, x=cluster_count_factor)) + geom_boxplot() + theme_bw() + facet_wrap(. ~ name, scales = "free_y")

dir.create(here(path_to_recip_hits, "coordinated_expression"))
#split up into three separate plots - Diagonal, percentile, top 1 - in that order of stringency. 
pdf(here(path_to_recip_hits, "coordinated_expression", "coordinated_expression_mean_Spearman.pdf"), height=5, width = 5)
ggplot(data = all_results_long %>% filter(name == "mean_diag_cor"), aes(y=value, x=cluster_count_factor)) + 
  geom_boxplot() + theme_bw() + xlab("Number of clusters") + ylab("Mean same gene Spearman correlation") +
  theme(axis.text.x = element_text(angle=270+45, hjust = 0))
dev.off()

pdf(here(path_to_recip_hits, "coordinated_expression", "coordinated_expression/coordinated_expression_mean_percentile_same_gene.pdf"), height=5, width = 5)
ggplot(data = all_results_long %>% filter(name == "mean_percentile_same_gene"), aes(y=value, x=cluster_count_factor)) + 
  geom_boxplot() + theme_bw() + xlab("Number of clusters") + ylab("Mean percentile of same gene correlation") +
  theme(axis.text.x = element_text(angle=270+45, hjust = 0))
dev.off()

pdf(here(path_to_recip_hits, "coordinated_expression", "coordinated_expression_top_1_proportion.pdf"), height=5, width = 5)
ggplot(data = all_results_long %>% filter(name == "top_1_proportion"), aes(y=value, x=cluster_count_factor)) + 
  geom_boxplot() + theme_bw() + xlab("Number of clusters") + ylab("Mean proportion of same gene top hits") +
  theme(axis.text.x = element_text(angle=270+45, hjust = 0))
dev.off()



