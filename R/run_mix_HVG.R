library(readr)
library(tibble)
library(zellkonverter)
library(Seurat)
library(mixhvg)
library(SingleCellExperiment)
library(scater)
library(here)

#hvg_subset <- readH5AD(here("data", "whole_mouse_brain", "processed", "merged_before_HVG.01pct.h5ad"), reader = "python", X_name="counts")
#loading 5% of the cells due to R memory limits
hvg_subset <- readH5AD(here("data/whole_mouse_brain/processed/merged_before_HVG.05pct.h5ad"), reader = "python", X_name="counts")

hvg_subset <- logNormCounts(hvg_subset)
seurat_obj <- as.Seurat(hvg_subset)
seurat_obj <- FindVariableFeaturesMix(seurat_obj, nfeatures = 3534)

head(VariableFeatures(seurat_obj))
seurat_obj


tibble(gene_id = VariableFeatures(seurat_obj))

tibble(gene_id = VariableFeatures(seurat_obj)) %>% write_csv(here("results", "HVG_lists", "mixhvg_FindVariableFeaturesMix_3534_5percent.csv"))

#write out, subset merged on high mem, run on narval
#run narval on original HVG subsets