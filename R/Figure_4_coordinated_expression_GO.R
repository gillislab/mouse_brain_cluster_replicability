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

#A key point is that reproduciblity of this is linked to the GO version used
GO.db
org.Mm.eg.db

#----------------------------
# GO children helper
#----------------------------
get_go_children <- function(go_id,
                            ontology = NULL,   # "BP", "MF", "CC" or NULL to auto-detect
                            recursive = TRUE,  # TRUE = all descendants; FALSE = direct children only
                            include_self = FALSE,
                            verbose = FALSE) {
  
  stopifnot(is.character(go_id), length(go_id) == 1)
  
  # Auto-detect ontology if not provided
  if (is.null(ontology)) {
    term_obj <- AnnotationDbi::mget(go_id, GO.db::GOTERM, ifnotfound = NA)[[1]]
    if (is(term_obj, "GOTerms")) {
      ontology <- term_obj@Ontology
    } else {
      if (verbose) message("Could not auto-detect ontology for ", go_id, " (not in GO.db snapshot?)")
      return(character(0))
    }
  }
  
  ont <- toupper(ontology)
  map <- switch(
    ont,
    BP = if (recursive) GO.db::GOBPOFFSPRING else GO.db::GOBPCHILDREN,
    MF = if (recursive) GO.db::GOMFOFFSPRING else GO.db::GOMFCHILDREN,
    CC = if (recursive) GO.db::GOCCOFFSPRING else GO.db::GOCCCHILDREN,
    stop("ontology must be one of 'BP', 'MF', 'CC' (or NULL to auto-detect).")
  )
  
  kids <- AnnotationDbi::mget(go_id, map, ifnotfound = NA)[[1]]
  if (length(kids) == 1 && is.na(kids)) kids <- character(0)
  
  kids <- unique(as.character(kids))
  if (include_self) kids <- unique(c(go_id, kids))
  kids
}


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
length(HVGs)
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
#----------------------------
# SynGO: expand to include children
#----------------------------
synGO_all_children <- character(0)
for (id in synGO_groups) {
  kids <- get_go_children(id, recursive = TRUE, include_self = FALSE, verbose = FALSE)
  if (length(kids)) synGO_all_children <- c(synGO_all_children, kids)
}
synGO_all_children <- unique(synGO_all_children)
length(synGO_all_children)

# final SynGO set: originals + descendants
synGO_groups <- unique(c(synGO_groups, synGO_all_children))
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
pdf(here(path_to_recip_hits, "coordinated_expression_top_1_proportion.pdf"), height=5, width = 5)
ggplot(data = to_plot, aes(y = avg_gene_exp_cor, x = proportion_recip, color = is_synGO)) + 
  geom_point(shape = 16, color = 'black') + theme_bw() + theme(legend.position = 'none') + theme(text = element_text(size = 18)) +
  geom_point(data = to_plot %>% filter(is_synGO), shape = 16, color = 'red') + scale_size_manual(values = c(0.5, 3)) +
  xlab("Proportion of reciprocal\nbest-hit genes") + ylab("Avg. gene expression correlation")
dev.off()
here(path_to_recip_hits, "coordinated_expression_top_1_proportion.pdf")

#note that this table may not match exactly as it's sensitive to GO versions
to_plot %>% select(GO_id, GO_name, is_synGO, num_genes = n_genes, average_correlation = avg_gene_exp_cor, prop_best = proportion_recip, gene_set = genes) %>% 
  write_csv(here(path_to_recip_hits, "coordinated_expression_GO_groups.csv"))

