library(dplyr)
library(ggplot2)
library(ggrepel)

studies = c('Ma_human', 'Franjic_human', 'Sepp_human', 'Ma_chimp', 'Ma_rhesus', 'Franjic_rhesus',
            'Ma_marmoset', 'Franjic_pig', 'Sepp_opossum', 'Tosches_lizard', 'Tosches_turtle',
            'Colquitt_zebrafinch', 'Colquitt_Bengalese_finch', 'Liao_frog', 'Lust_axolotl',
            'Anneser_zebrafish', 'Lamanna_lamprey', 'Davie_fruitfly', 'Jin_spider')


# load summary of mapping to Macosko/Zeng atlases
df1 = read.delim('cross_species_studies_mapping_summary_Macosko.csv', sep = ',')
df2 = read.delim('cross_species_studies_mapping_summary_Zeng.csv', sep = ',')

# get average AUROC per study x whole mouse brain atlas
t1 = df1 %>% group_by(target_study, species) %>% reframe(count = mean(AUROC))
t2 = df2 %>% group_by(target_study, species) %>% reframe(count = mean(AUROC))

t3 = data.frame(study = unlist(t1[,1]), species = unlist(t1[,2]), 
                Macosko = unlist(t1[,3]), Zeng = unlist(t2[,3]))
t3$avg_AUROC = 0.5*(t3$Zeng + t3$Macosko)
rownames(t3) = NULL


# get one AUROC per species --- averaged across studies
t4 = data.frame(t3 %>% group_by(species) %>% reframe(Zeng = mean(Zeng), Macosko = mean(Macosko),
                                     avg_AUROC = mean(avg_AUROC)))


# phylogenetic distance from TimeTree for plotting
spe_list = c('human', 'chimp', 'rhesus', 'marmoset', 'pig', 'opossum',
             'Bengalese_finch', 'zebrafinch', 'lizard', 'turtle', 'axolotl',
             'frog', 'zebrafish', 'lamprey', 'fruitfly', 'spider')

spe_ages = c(rep(87, 4), 94, 160, rep(319, 4), rep(352, 2),
            429, 563, rep(686, 2))

t4$phylo_time = spe_ages[match(t4$species, spe_list)]


# plot AUROC vs phylogenetic distance; using log-scaled (1-AUROC)
pdf('cross_species_auroc_by_age.pdf', width = 5, height = 5)
ggplot(t4, aes(x = phylo_time, y = 1-avg_AUROC, label = species)) + geom_point() +  
scale_y_log10(breaks = c(0.03, 0.1, 0.3), labels = c("0.97", "0.9", "0.7")) +
geom_text_repel() + theme_bw() + geom_abline(linetype = 'dashed', col = '#666666') + 
theme(text = element_text(size = 14)) + xlab('Phylo time (MYA)') + ylab('Avg AUROC') 
dev.off()