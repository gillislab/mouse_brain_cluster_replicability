library(SingleCellExperiment)
library(zellkonverter)
library(parallel)
library(magrittr)
library(dplyr)
library(here)
library(readr)
library(conflicted)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)

#IMPORTANT - this script needs switching for the two datasets by uncommenting code (this isn't ran often). That's at the 'switch here for Macosko/Zeng' lines

#create centroids for co-expression
#switch here for Macosko/Zeng
subset_path <- "/grid/gillis/data_norepl/lfrench/whole_mouse_brain/processed/macosko/subsets/"
#subset_path <- "/grid/gillis/data_norepl/lfrench/whole_mouse_brain/processed/zeng/subsets/"

average_log2_cpm <- NULL
#iterate parts
#sc.pp.normalize_total(adata_Zeng, target_sum=1e6, inplace=True) was used on the parts

#switch here for Macosko/Zeng
#for(part_file in list.files(subset_path, pattern = "AIT21.0.merged.with_multiome.cpm.part", full.names = TRUE)) {
for(part_file in list.files(subset_path, pattern = "Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.cpm.part_", full.names = TRUE)) {
  print(part_file)
  h5ad_part <- NULL
  gc()
  h5ad_part <- readH5AD(part_file, verbose=TRUE)
  #all_cell_types <-  h5ad_part$cl %>% unique()
  all_cell_types <-  h5ad_part$ClusterNm %>% unique()
  
  print("Done reading")
  
  log2_cpm_matrix <- log2(1+assay(h5ad_part, "X"))
  

  all_cell_types <- as.character(all_cell_types)
  average_expression_matrix <- matrix(NA, nrow = nrow(log2_cpm_matrix), ncol = length(all_cell_types),
                                      dimnames = list(rownames(log2_cpm_matrix), all_cell_types))
  print("Making dense")
  log2_cpm_matrix <- as.matrix(log2_cpm_matrix) #make it dense
  print("Done making dense")
  
  # Loop through each cell type and calculate average expression - could be faster with a matrix
  for (cell_type in all_cell_types) {
    if (grepl("000", cell_type)) {
      print(cell_type) #only works on Zeng
    }
    # Subset the SCE object for the current cell type
    #switch here for Macosko/Zeng
    #cell_type_subset <- log2_cpm_matrix[, h5ad_part$cl == cell_type]
    cell_type_subset <- log2_cpm_matrix[, h5ad_part$ClusterNm == cell_type]
    # Calculate average expression for each gene
    average_expression <- rowMeans(cell_type_subset)
    
    # Store the result in the matrix
    average_expression_matrix[, cell_type] <- average_expression
  }
  average_expression_matrix %<>% as.data.frame() %>% tibble()
  colnames(average_expression_matrix) <- all_cell_types
  average_expression_matrix %<>% mutate(gene_id = rownames(log2_cpm_matrix)) %>% select(gene_id, everything())
  
  average_log2_cpm <- bind_rows(average_log2_cpm, average_expression_matrix)
  print("Done one part, current size:")
  print(dim(average_log2_cpm))
}


#switch here for Macosko/Zeng
#write_csv(average_log2_cpm, path = paste0(subset_path, "AIT21.0.merged.with_multiome.log2_cpm.mean_centroids.csv.gz"))
write_csv(average_log2_cpm, path = paste0(subset_path, "Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.log2_cpm.mean_centroids.csv.gz"))

print("Done!")

