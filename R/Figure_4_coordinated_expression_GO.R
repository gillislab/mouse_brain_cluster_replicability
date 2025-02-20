library(ggplot2)
library(here)
library(readr)
library(readxl)
library(tidyr)
library(magrittr)
library(conflicted)
library(dplyr)
library(org.Mm.eg.db) #from bioconductor
library(GO.db)        #from bioconductor
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("rename", "dplyr")

path_to_recip_hits <- here("results","full_run_ZengAWS.1698256033/")

#provides a GO ID to NCBI gene ID mapping
goAnnots <- as.list(org.Mm.egGO2ALLEGS)
#convert to table in long format
goAnnots <- tibble(GO_id = names(goAnnots), NCBI_ID = goAnnots) %>% unnest(NCBI_ID)
#convert NCBI gene ID's to gene symbols
goAnnots <- left_join(goAnnots, org.Mm.egSYMBOL %>% as.data.frame() %>% rename(NCBI_ID = gene_id, gene_symbol = symbol))
nrow(goAnnots)
goAnnots <- inner_join(goAnnots, org.Mm.egENSEMBL %>% as.data.frame() %>% rename(NCBI_ID = gene_id) %>% distinct())
nrow(goAnnots)
#remove dupes
goAnnots %<>% distinct()

#filter for non HVGs
table_for_HVGs <- read_csv(here(path_to_recip_hits, "coordinated_expression_by_gene.csv"))
HVGs <- table_for_HVGs %>% filter(is_HVG) %>% pull(gene_id)
goAnnots %<>% filter(!ensembl_id %in% HVGs)

#add in aspect
goAnnots %<>% mutate(GO_aspect = Ontology(GO_id))
#add in GO group name
goAnnots %<>% mutate(GO_name = Term(GO_id))
#add in GO group size info
goAnnots %<>% inner_join(goAnnots %>% group_by(GO_id) %>% count() %>% rename(GO_group_size = n))

#look up a single gene
goAnnots %>% filter(gene_symbol == "Lig3")

#file from https://www.syngoportal.org/data/SynGO_bulk_download_release_20231201.zip
synGO_groups <- read_xlsx(here("data", "synGO" ,"release_20231201", "syngo_annotations.xlsx")) %>% pull(go_id) %>% unique()
length(synGO_groups)

#just BP
goAnnots %<>% filter(GO_aspect == "BP")

goAnnots %<>% inner_join(table_for_HVGs %>% select(ensembl_id = gene_id, cor_abs, is_same_gene_best))
goAnnots %>% pull(gene_symbol) %>% unique() %>% length()
to_plot <- goAnnots %>% group_by(GO_id, GO_name) %>% summarize(n_genes = n(), proportion_recip = sum(is_same_gene_best)/n(), 
                                                      avg_gene_exp_cor = mean(cor_abs), 
                                                      genes = paste0(gene_symbol, collapse = ','))
#were filtered to have between 15 and 150 genes
to_plot %<>% filter(n_genes <= 150, n_genes >= 15)
to_plot %<>% mutate(is_synGO = GO_id %in% synGO_groups)
to_plot %>% group_by(is_synGO) %>% count()

#Note - this will look different with different GO annotation database versions
pdf(here(path_to_recip_hits, "coordinated_expression", "coordinated_expression_top_1_proportion.pdf"), height=5, width = 5)
ggplot(data = to_plot, aes(y = avg_gene_exp_cor, x = proportion_recip, color = is_synGO)) + 
  geom_point() + theme_bw() + geom_point(data = to_plot %>% filter(is_synGO))
dev.off()

#note that this table may not match exactly as it's sensitive to GO versions
to_plot %>% select(GO_id, GO_name, is_synGO, num_genes = n_genes, average_correlation = avg_gene_exp_cor, prop_best = proportion_recip, gene_set = genes) %>% 
  write_csv(here(path_to_recip_hits, "coordinated_expression_GO_groups.csv"))

