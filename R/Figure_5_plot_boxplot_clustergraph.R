library(dplyr)
library(ggplot2)
library(igraph)

# load table with results of mapping cross-species clusters to
# Zeng/Macosko atlases
tab1 = read.delim('~/results/cross_species_studies/cross_species_mapping_results.tsv', sep = '\t')
tab1 <- tab1[!is.na(tab1$ref_celltype) & !is.na(tab1$target_celltype),]


# plot box plot across the three evolutionary conservation categories
plotdf = data.frame(AUROC = tab1$avg_AUROC,
                    fntype = tab1$evo_category)
plotdf <- plotdf[!is.na(plotdf$AUROC),]

pdf('plot-boxplot-evo-conservation.pdf', width = 5, height = 7)
ggplot(plotdf, aes(x = fntype, y = AUROC)) + geom_boxplot(outlier.shape = NA) +
theme_bw() + xlab('') + theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()



# ..... code snippet for plotting cluster graph ..... #
# get matrix of MetaNeighbor AUROC for Zeng/Macosko clusters vs cross-species clusters

# order of studies for matrix
studies_order = c('Ma_human', 'Franjic_human', 'Sepp_human', 'Ma_chimp', 'Ma_rhesus', 'Franjic_rhesus',
            'Ma_marmoset', 'Franjic_pig', 'Sepp_opossum', 'Tosches_lizard', 'Tosches_turtle',
            'Colquitt_zebrafinch', 'Colquitt_Bengalese_finch', 'Liao_frog', 'Lust_axolotl',
            'Anneser_zebrafish', 'Lamanna_lamprey', 'Davie_fruitfly', 'Jin_spider')

# list of cross-species datasets
dirs1 = c(rep('Tosches', 2), rep('Franjic', 3), rep('Sepp', 2), 'Lamanna', 'Lust',
         rep('Colquitt', 2), 'Anneser', rep('Ma', 4), 'Jin', 'Davie', 'Liao')
dirs2 = c('turtle', 'lizard', 'pig', 'rhesus', 'human', 'human', 'opossum', 'lamprey', 'axolotl',
         'zebrafinch', 'Bengalese_finch', 'zebrafish', 'human', 'chimp', 'rhesus', 'marmoset', 'spider',
         'fruitfly', 'frog')


# get MetaNeighbor scores for all datasets
pb = txtProgressBar(min = 0, max = length(dirs1), initial = 0)

auroc_mat = c()
for(id in 1:length(dirs1)){
    currdir = dirs1[id]
    currstudy = dirs2[id]
    
    aurocs1 = read.delim(paste0(currdir, '/', currstudy, '_pretrained_Macosko_full_MN_allvall.csv'),
                    sep = ',', check.names = F)
    aurocs2 = read.delim(paste0(currdir, '/', currstudy, '_pretrained_Zeng_full_MN_allvall.csv'),
                    sep = ',', check.names = F)

    temp = cbind(aurocs1, aurocs2)
    rownames(temp) <- paste0(currdir, '_', currstudy, '|', sub('.*\\|', '', rownames(temp)))
    
    auroc_mat = rbind(auroc_mat, temp)
    setTxtProgressBar(pb, id)
    
}


# subset to cross-species clusters mapped to Zeng/Macosko cell types with
# similar broad annotation
mat1 = tab1[!is.na(tab1$ref_level3) & (tab1$ref_level3==tab1$target_level3) & (tab1$is_recip==T) & (tab1$Zeng_celltype==tab1$Zeng_recip_celltype),]
mat1 <- mat1[!is.na(mat1$target_study),]

mat2 = mat1 %>% group_by(ref_celltype) %>% reframe(count = length(unique(target_study)), total = length(target_study)) %>% arrange(-count)
mat2$anno = tab1$ref_level2[match(mat2$ref_celltype, tab1$ref_celltype)]

# list of broad cell types in evolutionarily conserved clusters
ctypes = c('CGE InN', 'MGE InN', 'InN', 'ExN', 'Neural_progenitor',
           'Immune', 'Vascular', 'Astrocyte', 'Oligo', 'OPC')

mat2$anno_id = match(mat2$anno, ctypes)
mat3 <- arrange(mat2, anno_id)


# get table of cross-species clusters that map to one or both ends of 2,009 reciprocal hits
macs = mat3$ref_celltype
zeng = mat1$other_recip_celltype[match(mat3$ref_celltype, mat1$ref_celltype)]

list1 = c()
list2 = c()
df10 = c()

for(ii in 1:dim(mat3)[1]){
    list1 = c(list1, c(macs[ii], zeng[ii]))

    # order of cols
    order1 = order(match(mat1$target_study[mat1$ref_celltype==macs[ii]], studies_order))
    cls1 = mat1$target_celltype[mat1$ref_celltype==macs[ii]]
    list2 = c(list2, cls1[order1])

    # reshape mat3 into long table
    temp2 = do.call("rbind", replicate(length(cls1), mat3[ii,], simplify = FALSE)) 
    temp2$target_celltype = cls1[order1]
    df10 = rbind(df10, temp2)
    
}


# get adjacency matrix for cluster graph
labels1 = c(macs, list2)
adj = matrix(0, nrow = length(labels1), ncol = length(labels1))

for(ii in 1:dim(df10)[1]){
    i1 = match(df10$ref_celltype[ii], labels1)
    i2 = match(df10$target_celltype[ii], labels1)
    
    adj[i1, i2] = 1
    adj[i2, i1] = 1
}
rownames(adj) = labels1
colnames(adj) = labels1

g1 <- graph_from_adjacency_matrix(adj, mode = 'undirected')


# plot cluster graph
pdf('cross_species_cluster-graph.pdf', width = 20, height = 20)
plot(g1, vertex.label.color = 'white', vertex.size = 2, vertex.label.dist = 0.5)
dev.off()
