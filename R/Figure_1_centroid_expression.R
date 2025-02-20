library(gplots)
library(ggplot2)
library(data.table)
library(matrixStats)
library(magrittr)
library(dplyr)
library(tidyr)
library(here)
library(viridis)
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

single_cell_centroids <- read_csv(paste0(subset_path_cell, "AIT21.0.merged.with_multiome.log2_cpm.mean_centroids.csv.gz"))
single_cell_centroids[1:5,1:5]
dim(single_cell_centroids)

#############
#convert to a matrix
single_cell_centroids_matrix <- as.matrix(single_cell_centroids %>% select(-gene_id))
rownames(single_cell_centroids_matrix) <- single_cell_centroids %>% pull(gene_id)
dim(single_cell_centroids_matrix)
single_cell_centroids_matrix_full <- single_cell_centroids_matrix


single_nuc_centroids <- read_csv(paste0(subset_path_nuc, "Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.log2_cpm.mean_centroids.csv.gz"))

single_nuc_centroids_matrix <- as.matrix(single_nuc_centroids %>% select(-gene_id))
rownames(single_nuc_centroids_matrix) <- single_nuc_centroids %>% pull(gene_id)
intersect(rownames(single_nuc_centroids_matrix), single_cell_centroids %>% pull(gene_id)) %>% length() #21649
#single_nuc_centroids_matrix[1:5,1:5]
#dim(single_nuc_centroids_matrix)
single_nuc_centroids_matrix_full <- single_nuc_centroids_matrix #full version before 2009 set

single_nuc_centroids_matrix[1:5,1:5]
single_cell_centroids_matrix[1:5,1:5]

common_genes <- intersect(rownames(single_nuc_centroids_matrix), rownames(single_cell_centroids_matrix))
length(common_genes)

single_nuc_centroids_matrix_full <- single_nuc_centroids_matrix_full[common_genes, ]
single_cell_centroids_matrix_full <- single_cell_centroids_matrix_full[common_genes, ]


gene_means_cell <- matrixStats::rowMeans2(single_cell_centroids_matrix_full)
gene_means_nuc <- matrixStats::rowMeans2(single_nuc_centroids_matrix_full)

gene_means_cell <- tibble(gene_id = rownames(single_cell_centroids_matrix_full), mean_exp = gene_means_cell)
gene_means_nuc <- tibble(gene_id = rownames(single_nuc_centroids_matrix_full), mean_exp = gene_means_nuc)
cor.test(gene_means_cell$mean_exp, gene_means_nuc$mean_exp, m='s')
cor.test(gene_means_cell$mean_exp, gene_means_nuc$mean_exp, m='p')

gene_means_joined <- inner_join(gene_means_cell, gene_means_nuc, by="gene_id") %>% rename(mean_exp_single_cell = mean_exp.x, mean_exp_single_nuc = mean_exp.y)
correlation <- cor(gene_means_joined$mean_exp_single_cell, gene_means_joined$mean_exp_single_nuc, m='s')
gene_means_joined %<>% mutate(distance_from_identity = abs(mean_exp_single_cell - mean_exp_single_nuc)/sqrt(2))
gene_means_joined %<>% mutate(directioned_distance = distance_from_identity * sign(mean_exp_single_cell - mean_exp_single_nuc))
gene_means_joined %>% arrange(-directioned_distance)
gene_means_joined %>% arrange(directioned_distance)

pdf(here("results", "average_expression_centroids.pdf"), height=6, width = 6)
ggplot(data=gene_means_joined, aes(x = mean_exp_single_cell, y = mean_exp_single_nuc)) + geom_point(alpha  = 0.2) + 
  theme_bw() + xlab("Average expression in cells (log(CPM))") +
  ylab("Average expression in nuclei (log(CPM))") + coord_fixed() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  annotate("text",
           x = 0.05, 
           y = 12, 
           label = paste0("Spearman r: ", correlation %>% format(digits = 2)),
           hjust = 0, 
           vjust = 1
  ) 
dev.off()

#for testing
#single_nuc_centroids_matrix_full <- single_nuc_centroids_matrix_full[1:300, 1:200]
#single_cell_centroids_matrix_full<- single_cell_centroids_matrix_full[1:300, 1:300]

#genome wide spearman correlation heatmap - takes awhile
mat1 <- cor(cbind(single_nuc_centroids_matrix_full, single_cell_centroids_matrix_full), m='s')
dim(mat1)
mat1[1:5,1:5]

hm <- heatmap.2(mat1, col = viridis(100), labCol = "", labRow = "",
                density.info = "none", trace = "none", scale = "none")

png(here("results", "centroids_correlation_heatmap.png"), width=14000, height = 14000)
row_order <- hm$rowInd
col_order <- hm$colInd

# Replot using stored ordering (this skips reclustering)
heatmap.2(mat1[row_order, col_order], col = viridis(100),
          labCol = "", labRow = "",
          density.info = "none", trace = "none", scale = "none",
          Rowv = FALSE, Colv = FALSE)
dev.off()

# Get the row order based on the dendrogram
row_order = order.dendrogram(hm$rowDendrogram)

# Get the rownames in the order they appear in the heatmap
ordered_rownames = rownames(mat1)[row_order]

newdf <- tibble(cluster = ordered_rownames) %>% 
  mutate(x = -1*row_number(), study = if_else(cluster %in% colnames(single_nuc_centroids_matrix_full), "Macosko", "Zeng"))
newdf %<>% mutate(color = recode(study, Macosko = '#E21E25', Zeng = '#4450A2'))
newdf %<>% mutate(study_type = recode(study, Macosko = "Nuclei", Zeng = "Cells"))
newdf %>% select(study_type, study) %>% distinct()
pdf(here("results", "average_expression_centroids_color_bar.pdf"), height=30, width = 4)
ggplot(newdf, aes(x = 1, y = x, fill = study_type)) +
  geom_tile() +
  scale_fill_manual("Dataset", values = setNames(newdf %>% pull(color), newdf %>% pull(study_type)), na.value="white") +
  theme_minimal() + ylab("") + xlab("") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())
dev.off()
