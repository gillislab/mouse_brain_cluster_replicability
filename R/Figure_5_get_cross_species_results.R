library(dplyr)

# load summary of top hits from Macosko/Zeng pretained models
df1 = read.delim('~/results/cross_species_studies/cross_species_studies_mapping_summary_Macosko.csv', sep = ',')
df2 = read.delim('~/results/cross_species_studies/cross_species_studies_mapping_summary_Zeng.csv', sep = ',')

df1$ref = 'Macosko'
df2$ref = 'Zeng'

# load table with broad cell type annotations for
# 5k clusters in Zeng and Macosko atlases
newdf = read.delim('~/results/cross_species_studies/Zeng_Macosko_broad_celltype_annotation.tsv', sep = '\t')

# load table of broad cell type annotations across sampled species
tab1 = read.delim('~/results/cross_species_studies/cross_species_cluster_annotations.csv', sep = ',')


# ... add broad cell type annotations to Zeng/Macosko top hits summaries ... #
df1$target_level1 = tab1$level1_anno[match(df1$target_celltype, tab1$celltype)]
df1$target_level2 = tab1$level2_anno[match(df1$target_celltype, tab1$celltype)]

# making coarse resolution annotations for Macosko
df1$target_level3 = df1$target_level2
df1$target_level3[df1$target_level1=='Neuron'] = 'ExN'
df1$target_level3[grep('InN', df1$target_level2)] = 'InN'
df1$target_level3[grep('Astrocyte|Ependymal', df1$target_level2)] = 'AstroEpendymal'


df1$ref_level1 = newdf$level1_anno[match(df1$ref_celltype, newdf$cluster)]
df1$ref_level2 = newdf$level2_anno[match(df1$ref_celltype, newdf$cluster)]

df1$ref_level3 = df1$ref_level2
df1$ref_level3[df1$ref_level1=='Neuron'] = 'ExN'
df1$ref_level3[grep('InN', df1$ref_level2)] = 'InN'
df1$ref_level3[grep('Astrocyte|Ependymal', df1$ref_level2)] = 'AstroEpendymal'


df2$target_level1 = tab1$level1_anno[match(df2$target_celltype, tab1$celltype)]
df2$target_level2 = tab1$level2_anno[match(df2$target_celltype, tab1$celltype)]

# making coarse resolution annotations for Zeng
df2$target_level3 = df2$target_level2
df2$target_level3[df2$target_level1=='Neuron'] = 'ExN'
df2$target_level3[grep('InN', df2$target_level2)] = 'InN'
df2$target_level3[grep('Astrocyte|Ependymal', df2$target_level2)] = 'AstroEpendymal'


df2$ref_level1 = newdf$level1_anno[match(df2$ref_celltype, newdf$cluster)]
df2$ref_level2 = newdf$level2_anno[match(df2$ref_celltype, newdf$cluster)]

df2$ref_level3 = df2$ref_level2
df2$ref_level3[df2$ref_level1=='Neuron'] = 'ExN'
df2$ref_level3[grep('InN', df2$ref_level2)] = 'InN'
df2$ref_level3[grep('Astrocyte|Ependymal', df2$ref_level2)] = 'AstroEpendymal'


# exclude these frog clusters since they're not sampled by the atlases
# this step reduces 1,204 cross-species clusters to 1,192 clusters 
exc_frog_types = c('Endocrine progenitor cell_chgb high', 'Endocrine progenitor cell_st6galnac2 high',
       'Gonadotroph cell_cga high', 'Gonadotroph cell_cyp21a2.1 high', 'Growth hormone cell_Rp high',
       'Growth hormone cell_hnrnpk high', 'Growth hormone cell_rasd1 high', 'Growth hormone progenitor cell',
       'Melanotrope cell', 'Rp_high', 'Thyrotroph cell', 'Prolactin progenitor cell')

newdf1 <- df1[- which(df1$target_celltype %in% paste0('Liao_frog|', exc_frog_types)),]
newdf2 <- df2[- which(df2$target_celltype %in% paste0('Liao_frog|', exc_frog_types)),]



## evolutionarily conserved cell types? ##
# how many matched to reciprocal and correctly predicted?
mat1 = newdf1[!is.na(newdf1$ref_level3) & (newdf1$ref_level3==newdf1$target_level3) & (newdf1$is_recip==T),]

mat2 = newdf2[!is.na(newdf2$ref_level3) & (newdf2$ref_level3==newdf2$target_level3) & (newdf2$is_recip==T),]


# load table with Zeng and Macosko clusters and their reciprocal hit (or NA)
# in the other study
anno = read.delim('~/results/cross_species_studies/mouse_brain_cluster_map_summary.csv', sep = ',')



# do cross-species clusters match to both sides of 2,009 pairs?
newdf1$other_celltype = newdf2$ref_celltype[match(newdf1$target_celltype, newdf2$target_celltype)]
newdf1$other_recip = anno$is_reciprocal[match(newdf1$other_celltype, anno$cluster)]
newdf1$other_recip_celltype = anno$recip_cluster[match(newdf1$ref_celltype, anno$cluster)]

newdf1$other_AUROC = newdf2$AUROC[match(newdf1$target_celltype, newdf2$target_celltype)]
newdf1$avg_AUROC = 0.5*(newdf1$AUROC + newdf1$other_AUROC)


newdf2$other_celltype = newdf1$ref_celltype[match(newdf2$target_celltype, newdf1$target_celltype)]
newdf2$other_recip = anno$is_reciprocal[match(newdf2$other_celltype, anno$cluster)]
newdf2$other_recip_celltype = anno$recip_cluster[match(newdf2$ref_celltype, anno$cluster)]

newdf2$other_AUROC = newdf1$AUROC[match(newdf2$target_celltype, newdf1$target_celltype)]
newdf2$avg_AUROC = 0.5*(newdf2$AUROC + newdf2$other_AUROC)


# stratify cross-species clusters based on whether they map to both ends of 2,009 pairs
# one end of 2,009 pairs, or none of the 2,009 reciprocal clusters

# clusters hitting both ends of 2,009 recip pairs
c1 = newdf1$ref_celltype[(!is.na(newdf1$ref_level3) & (newdf1$ref_level3==newdf1$target_level3) & (newdf1$is_recip==T) & (newdf1$other_celltype==newdf1$other_recip_celltype))]
c2 = newdf2$ref_celltype[(!is.na(newdf2$ref_level3) & (newdf2$ref_level3==newdf2$target_level3) & (newdf2$is_recip==T) & (newdf2$other_celltype==newdf2$other_recip_celltype))]
c1 <- c1[!is.na(c1)]
c2 <- c2[!is.na(c2)]

# clusters hitting either of 2,009 recip pair
c3 = newdf1$ref_celltype[(!is.na(newdf1$ref_level3) & (newdf1$ref_level3==newdf1$target_level3) & (newdf1$is_recip==T))]
c4 = newdf2$ref_celltype[(!is.na(newdf2$ref_level3) & (newdf2$ref_level3==newdf2$target_level3) & (newdf2$is_recip==T))]
c3 <- c3[!is.na(c3)]
c4 <- c4[!is.na(c4)]

# unique ref cell types that only match to Macosko or Zeng recip
c5 = setdiff(unique(union(c3, c4)), unique(union(c1, c2)))


# create data.frame to plot the distribution of AUROCs for
# cross-species clusters in the three categories
combodf = newdf1
combodf$evo_category = 'diverged'

combodf$evo_category[combodf$ref_celltype %in% unique(c1)] = 'conserved_both'
combodf$evo_category[combodf$ref_celltype %in% c5 | combodf$other_celltype %in% c5] = 'conserved_one'
colnames(combodf)[15:18] = c('Zeng_celltype', 'Zeng_recip', 'Zeng_recip_celltype', 'Zeng_AUROC')


# save combodf for plotting evo conserved cell types
cols1 = c('target_study', 'species', 'target_celltype', 'target_size', 'target_level1', 'target_level2', 'target_level3',
          'ref_celltype', 'is_recip', 'Zeng_recip_celltype', 'AUROC', 'ref_level1', 'ref_level2', 'ref_level3',
          'Zeng_celltype', 'Zeng_recip', 'Zeng_AUROC', 'avg_AUROC', 'evo_category')

write.table(combodf[,cols1], file = '~/results/cross_species_studies/cross_species_mapping_results.tsv', sep = '\t', row.names = F,
           col.names = T, quote = F)
