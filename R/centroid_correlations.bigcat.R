#analysis code is in exploratory/examine_MERFISH_calls_for_recips
library(SingleCellExperiment)
library(zellkonverter)
library(parallel)
#remotes::install_github("AllenInstitute/scrattch.bigcat")
library(scrattch.bigcat)
library(MetaMarkers)
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
conflicts_prefer(scrattch.bigcat::cpm)

#load Macosko medians, convert to matrix, test get_cl_means?

#load up spatial data - convert to log transformed CPM and sparse matrix for R - zelkconverter may need deletion of .cache/R/basilisk/ before running again
Zeng_MERFISH <- readH5AD(here("data", "whole_mouse_brain",  "zeng", "MERFISH-C57BL6J-638850", "20230630", "C57BL6J-638850-raw.h5ad"), reader = "python", X_name="counts")

#write out base cell data for downstream scripts
h5ad_base <- colData(Zeng_MERFISH) %>% as.data.frame() %>% tibble()
h5ad_base %<>% mutate(cell_label = colnames(Zeng_MERFISH))
h5ad_base %<>% select(cell_label, everything())
h5ad_base %>% write_csv(here("results", "MERFISH500", "H5ad_calls.csv.gz"))

z <- counts(Zeng_MERFISH)
z[1:5,1:5]
SingleCellExperiment::cpm(Zeng_MERFISH) <- log2(1+cpm(z))

Sys.setenv("VROOM_CONNECTION_SIZE"= 21000000)
#from the google drive file linked from braincelldata.org
single_nuc_centroids <- read_csv(here("data", "whole_mouse_brain", "macosko", "from_google_drive", "Single_Nuc_Cluster_Avg_Expression.csv.gz"))
dim(single_nuc_centroids)

single_nuc_centroids %<>% rename(cluster_id = `...1`)

#need a:
#A gene (rows) x samples (columns) sparse matrix
test <- single_nuc_centroids
test_matrix <- test %>% select(-cluster_id) %>% as.matrix()
test_matrix[1:5,1:5]
rownames(test_matrix) <- test %>% pull(cluster_id)
test_matrix <- t(test_matrix)
test_matrix[1:5,1:5]
#####
#####
single_nuc_genes <- tibble(rowname = rownames(test_matrix))
single_nuc_genes %<>% mutate(gene_symbol = gsub("=.*","", rowname))
rowData(Zeng_MERFISH) %>% as.data.frame() %>% tibble() %>% pull(gene_symbol)
MERFISH_genes <- rowData(Zeng_MERFISH) %>% as.data.frame() %>% tibble() %>% pull(gene_symbol)
single_nuc_genes %<>% filter(gene_symbol %in% MERFISH_genes)
test_matrix <- test_matrix[single_nuc_genes %>% pull(rowname),]
dim(test_matrix)
rownames(test_matrix) <- gsub("=.*","", rownames(test_matrix))
test_matrix[1:5,1:5]

#call cpm from scrattch.bigcat
#the Macosko data is average counts which is further processed, lets use the same - or just use Spearman?
test_matrix <- log2(1+cpm(test_matrix))

#filter for spatial genes and deal with dupes
####
Zeng_cpm_matrix <- SingleCellExperiment::cpm( Zeng_MERFISH)
rownames(Zeng_cpm_matrix) <- MERFISH_genes
Zeng_cpm_matrix[1:5,1:5]

################################
################################

train.cl.dat <- test_matrix
test.dat <- Zeng_cpm_matrix[rownames(test_matrix),]
dim(train.cl.dat)
dim(test.dat)
saveRDS(test.dat, file = paste0(tempdir(), "/test.dat.rds"))
saveRDS(train.cl.dat, file = paste0(tempdir(), "/train.cl.dat.rds"))
tempdir()
rm(list = ls()) #clear workspace - should free up memory
gc()
gc()

#####################
#####################
#output folders
dir.create(here("results", "MERFISH500"))
dir.create(here("results", "MERFISH500", "Macosko_predictions"))

#####################
#####################
#reload from here
markers.perc = 0.8
iter = 100
verbose = TRUE
mc.cores=50

test.dat <- readRDS(file = paste0(tempdir(), "/test.dat.rds"))
train.cl.dat <- readRDS(file = paste0(tempdir(), "/train.cl.dat.rds"))
markers <- rownames(test.dat)
#markers <- head(markers, n=10)
test.dat <- test.dat[markers,]
#test.dat <- test.dat[,1:60000] #for testing

#cut it into 50 parts, run each and write out
number_of_splits <- 220
cells_to_split <- colnames(test.dat)
split_list <- split(cells_to_split, rep(1:number_of_splits, length.out = length(cells_to_split)))

#most code from scrattch.bigcat/R/annotate.R
#https://github.com/AllenInstitute/scrattch.bigcat/blob/bbfe5277de63c65a5342e43f3e1eeb1015755fc7/R/annotate.R

for(split_index in 1:length(split_list)) {
  print(paste0("Running index:", split_index))
  test.dat <- readRDS(file = paste0(tempdir(), "/test.dat.rds"))
  test.dat <- test.dat[, split_list[[split_index]]]
  gc()
  
  
  # Perform mapping iter times.
  #system.time({
  map.result <- mclapply(1:iter, 
                         function(i){
                           if(verbose) {
                             cat("\r", paste0("Running iteration ",i," of ",iter,".        "))
                             flush.console()
                           }
                           
                           #set the seed so every cell gets the same shuffles of the marker genes
                           set.seed(i)
                           
                           tmp.markers <- sample(markers, round(length(markers) * markers.perc))
                           test.cl.cor <- cor(as.matrix(test.dat[tmp.markers,]), train.cl.dat[tmp.markers,])
                           test.cl.cor[is.na(test.cl.cor)] <- 0
                           max.cl.cor <- apply(test.cl.cor, 1, which.max)
                           pred.cl <- colnames(test.cl.cor)[max.cl.cor]
                           pred.cl <- setNames(pred.cl, row.names(test.cl.cor))
                           pred.score <- apply(test.cl.cor, 1, max)
                           # if (is.factor(train.cl)) {
                           #   pred.cl <- setNames(factor(pred.cl, levels = levels(train.cl)), 
                           #                       names(pred.cl))
                           # }
                           pred.df <- data.frame(pred.cl = pred.cl, pred.score = pred.score)
                         }, mc.cores=mc.cores)
  #})
  
  # Extract predicted cluster assignments from each iteration
  map.cl <- sapply(map.result, 
                   function(x) {
                     x$pred.cl
                   })


    
  # Compute fraction of times each sample mapped to each cluster
  row.names(map.cl) <- colnames(test.dat)
  map <- as.data.frame(as.table(as.matrix(map.cl)))
  map.table <- table(map$Var1, map$Freq)
  map.freq <- unclass(map.table)
  
  # Find the most frequently mapped cluster for each sample
  max.freq <- apply(map.freq, 1, which.max)
  pred.cl <- colnames(map.freq)[max.freq]
  pred.cl <- setNames(pred.cl, row.names(map.freq))
  
  # Gather results
  map.df <- data.frame(pred.cl = pred.cl, 
                       prob = matrixStats::rowMaxs(map.freq) / iter)
  
  # output results
  out_list <- list(map.df = map.df,
                   map.freq = map.freq)
  
  #compute the score
  map.full <- lapply(map.result, as_tibble, rownames="cell_id")
  map.full <- bind_rows(map.full)
  
  output_mapping <- as_tibble(out_list$map.df, rownames="cell_id")
  joined_scores <- inner_join(map.full, output_mapping)
  joined_scores %<>% group_by(cell_id, pred.cl) %>% summarize(average_correlation_score = mean(pred.score), prob=mean(prob))

  write_csv(joined_scores, here("results", "MERFISH500", "Macosko_predictions", paste0("joined_scores.part_",split_index,".csv") ))
}
