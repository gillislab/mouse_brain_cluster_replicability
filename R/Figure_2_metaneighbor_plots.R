library(circlize)
library(plotROC)
library(viridis)
library(zoo)
library(dendextend)
library(tidyr)
library(ggpointdensity)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(readxl)
library(here)
library(magrittr)
library(MetaNeighbor)
library(reshape2)
library(conflicted)
conflicts_prefer(dplyr::last)
conflict_prefer("first", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("union", "base")
conflict_prefer("setdiff", "base")
Sys.setenv("VROOM_CONNECTION_SIZE"=1000000)

base_folder <- here("results/full_run_ZengAWS.1698256033/") #updated Allen data - 2009 recips

all_AUROCs <- read_csv(paste0(base_folder, "aurocs_full.csv.gz"))
all_AUROCs %<>% rename(target = "...1")
dim(all_AUROCs)
all_AUROCs_matrix <- as.matrix(all_AUROCs %>% select(-target))
rownames(all_AUROCs_matrix) <- all_AUROCs %>% pull(target)

top_hits <- read_csv(paste0(base_folder, "top_hits.0.95.csv"))
recip_types <- union(top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_1`),
                     top_hits %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(`Study_ID|Celltype_2`))

one_vs_one_AUROCs <- read_csv(paste0(base_folder, "aurocs_1v1.csv.gz"), guess_max = 200000000)
one_vs_one_AUROCs %<>% rename(target = "...1")
one_vs_one_AUROCs_melted <- melt(one_vs_one_AUROCs, na.rm = T, value.name= "auroc", variable.name = "reference") %>% tibble()
one_vs_one_AUROCs_melted %<>% mutate(reference = as.character(reference), target = as.character(target))
one_vs_one_AUROCs_melted %<>% mutate(reference_study = getStudyId(reference), target_study = getStudyId(target))

#check that the reference and target are not flipped - should be four for each reference
one_vs_one_AUROCs_melted %>% group_by(reference) %>% count() 

#convert AUROC's to p-values with past function
get_p_from_AUC <- function(n.x, n.y, AUC, alternative = "two.sided", correct = T) {
  n.x <- as.double(n.x)
  n.y <- as.double(n.y)
  
  STATISTIC = AUC*(n.x*n.y)
  z <- STATISTIC - n.x * n.y/2
  SIGMA <- sqrt((n.x * n.y/12) * ((n.x + n.y + 1) ))
  
  if (correct) {
    CORRECTION <- switch(alternative, two.sided = sign(z) * 0.5, greater = 0.5, less = -0.5)
  } else {
    CORRECTION <- 0
  }
  z <- (z - CORRECTION)/SIGMA
  PVAL <- switch(alternative, less = pnorm(z), greater = pnorm(z, lower.tail = FALSE), two.sided = 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE)))
  PVAL
}

merged_obs <- read_csv(paste0(base_folder, "merged.obs.csv.zip"))

# Get cell type sizes
cell_type_sizes <- merged_obs %>%
  group_by(cell.type = paste0(study_id, "|", cell.type)) %>%
  count()


one_vs_one_AUROCs_melted %<>% inner_join(cell_type_sizes %>% select(target = cell.type, target_size = n))
one_vs_one_AUROCs_melted %<>% mutate(within_study = reference_study == target_study)

one_vs_one_p_values <- one_vs_one_AUROCs_melted %>% group_by(reference, within_study) %>% 
  summarize(auroc = first(auroc), size1 = first(target_size), size2 = last(target_size)) %>% rowwise() %>% 
  mutate(p_value = get_p_from_AUC(size1, size2, auroc))

one_vs_one_p_values <- inner_join(one_vs_one_p_values, cell_type_sizes %>% select(reference = cell.type, reference_size = n))

one_vs_one_p_values %<>% mutate(reference_study = getStudyId(as.character(reference)))

one_vs_one_p_values %<>% mutate(comparison_type = recode(as.character(within_study), "TRUE" = "within study", "FALSE" = "across study"))
one_vs_one_p_values %<>% ungroup()
one_vs_one_p_values %<>% mutate(comparison_type = factor(comparison_type, levels = rev(sort(unique(comparison_type)))))
one_vs_one_p_values %<>% mutate(neg_log_p = -1*log10(p_value))
one_vs_one_p_values %>% arrange(-neg_log_p)
one_vs_one_p_values %<>% mutate(neg_log_p = ifelse(neg_log_p == Inf, max(neg_log_p[neg_log_p != Inf], na.rm = TRUE), neg_log_p))

one_vs_one_p_values %<>% mutate(`Reference Dataset` = if_else(reference_study == "Zeng", "Cells", "Nuclei"))

#supplement Figure 1
pdf(paste0(base_folder, "log_p_best_vs_next.pdf"), height=5, width = 6.5)
ggplot(one_vs_one_p_values, aes(y=neg_log_p, x = log10(reference_size), color = `Reference Dataset`, linetype = comparison_type)) + theme_bw() + 
  ylab("-log10(best vs next p-value)") +  xlab("log10(cluster size)") + 
  scale_color_manual(values = c("#4450A2", "#E21E25")) +
  geom_smooth() 
dev.off()

###########
###########
#Next is the density plots
###########
all_vs_all_AUROCs_melted <- melt(all_AUROCs_matrix, na.rm = T, value.name= "auroc", variable.name = "reference") %>% tibble()
all_vs_all_AUROCs_melted %<>% rename(target = Var1, reference = Var2, auroc_one_v_all = auroc)
all_vs_all_AUROCs_melted %<>% mutate(target = as.character(target), reference = as.character(reference))

for_plot <- inner_join(one_vs_one_AUROCs_melted %>% select(target, reference, auroc), all_vs_all_AUROCs_melted)
for_plot %<>% mutate(reference_study = getStudyId(as.character(reference)))
for_plot %<>% mutate(target_study = getStudyId(as.character(target)))


for_plot %<>% filter(reference_study != target_study)

#use the one with the highest auroc_one_v_all
for_plot %<>% group_by(reference) %>% arrange(desc(auroc_one_v_all)) %>% dplyr::slice(1)
for_plot %>% pull(auroc) %>% median()
for_plot %>% pull(auroc) %>% mean()

for_plot %>% group_by(reference_study) %>% summarize(mean_all = mean(auroc_one_v_all), mean_best_vs_next = mean(auroc), n=n())

for_plot %<>% mutate(log_inverse = -1*log10(1-auroc_one_v_all))
for_plot %<>% mutate(inverse = 1-auroc_one_v_all)

cor.test(for_plot$auroc, for_plot$auroc_one_v_all, m='s')
convert_ticks <- function(x) {
  1-(1/10**x)
}

pdf(paste0(base_folder, "Point_density_two_AUROCs_Zeng_reference.pdf"), width=6, height = 5)
(Zeng_plot <- ggplot(data = for_plot %>% filter(reference_study == "Zeng"), 
                     aes(x = log_inverse, y = auroc)) + 
    geom_pointdensity(adjust  = 0.1, size = 1) + 
    theme_bw() + 
    xlab("One vs all AUROC (log scale)") + 
    ylab("Best vs next AUROC") + labs(color = "Neighbors") +
    scale_color_viridis_c(limits = c(1,330)) + 
    scale_x_continuous(labels = convert_ticks) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))) #+
    #geom_hline(yintercept = 0.6) + geom_vline(xintercept = 2)) 
dev.off()

x_breaks <- ggplot_build(Zeng_plot)$layout$panel_params[[1]]$x$breaks

pdf(paste0(base_folder, "Point_density_two_AUROCs_Macosko_reference.pdf"), width=6, height = 5)
ggplot(data = for_plot %>% filter(reference_study == "Macosko"), 
       aes(x = log_inverse, y = auroc)) + 
  geom_pointdensity(adjust  = 0.1, size = 1) + 
  theme_bw() + 
  xlab("One vs all AUROC (log scale)") + 
  ylab("Best vs next AUROC")  + labs(color = "Neighbors") +
  scale_color_viridis_c(limits = c(1,330)) + 
  scale_x_continuous(labels = convert_ticks, breaks = x_breaks) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
  #geom_hline(yintercept = 0.6) + geom_vline(xintercept = 2)
dev.off()

##########
#Histograms
##########

#for AUROC plot, reuses the above
for_plot %>% pull(auroc) %>% mean() 
for_plot %<>% mutate(is_ref_recip = reference %in% recip_types)

for_plot %<>% mutate(is_ref_recip_string = if_else(is_ref_recip, "Reciprocal","Nonreciprocal"))

for_plot %>% group_by(is_ref_recip, reference_study) %>% count()
for_plot %>% group_by(is_ref_recip, reference_study) %>% summarize(n=n(), mean=mean(auroc))

for_plot_means <- for_plot %>% 
  group_by(is_ref_recip_string, reference_study) %>% 
  summarise(mean_auroc = mean(auroc))

#add bar - mean, make transparent, write out as pdf
pdf(paste0(base_folder, "Best_vs_next_AUROCs_Macosko_reference.pdf"), width=6.5, height = 5)
original_color <- "#E21E25" 
palette <- colorRampPalette(c(original_color, "#FFFFFF"))
palette_dark <- colorRampPalette(c(original_color, "#000000"))
ggplot(for_plot %>% filter(reference_study == "Macosko"), aes(x=auroc, fill = is_ref_recip_string)) + 
  geom_histogram(data= for_plot %>% filter(reference_study == "Macosko", is_ref_recip), alpha = 0.7) + 
  geom_histogram(data= for_plot %>% filter(reference_study == "Macosko", !is_ref_recip), alpha = 0.7) + 
  theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + xlab("Best vs next AUROC") +
  scale_fill_manual(values = c(palette_dark(4)[2], palette(4)[2])) + ylab("Frequency") + labs(fill = "Reference status") +
  geom_vline(data = for_plot_means %>% filter(reference_study == "Macosko"), aes(xintercept = mean_auroc, color = is_ref_recip_string), linetype = "dashed", linewidth = 1, color = c(palette_dark(4)[2], palette(4)[2]))
dev.off()

pdf(paste0(base_folder, "Best_vs_next_AUROCs_Zeng_reference.pdf"), width=6.5, height = 5)
original_color <- "#4450A2" 
palette <- colorRampPalette(c(original_color, "#FFFFFF"))
palette_dark <- colorRampPalette(c(original_color, "#000000"))
ggplot(for_plot %>% filter(reference_study == "Zeng"), aes(x=auroc, fill = is_ref_recip_string)) + 
  geom_histogram(data= for_plot %>% filter(reference_study == "Zeng", is_ref_recip), alpha = 0.7) + 
  geom_histogram(data= for_plot %>% filter(reference_study == "Zeng", !is_ref_recip), alpha = 0.7) + 
  theme_bw() + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + xlab("Best vs next AUROC") +
  scale_fill_manual(values = c(palette_dark(4)[2], palette(4)[2])) + ylab("Frequency") + labs(fill = "Reference status") +
  geom_vline(data = for_plot_means %>% filter(reference_study == "Zeng"), aes(xintercept = mean_auroc, color = is_ref_recip_string), linetype = "dashed", size = 1, color = c(palette_dark(4)[2], palette(4)[2]))
dev.off()

#Barplot of sizes
cell_type_sizes %<>% mutate(study_id = getStudyId(cell.type))

for_plot <- cell_type_sizes %>% group_by(study_id) %>% arrange(-n) %>% 
  mutate(cumulative_cell_fraction = cumsum(n)/sum(n), cumulative_cluster_fraction  = row_number()/n())
for_plot %<>% bind_rows(tibble(study_id = "Zeng", cumulative_cell_fraction = 0, cumulative_cluster_fraction =0))
for_plot %<>% bind_rows(tibble(study_id = "Macosko", cumulative_cell_fraction = 0, cumulative_cluster_fraction =0))
for_plot %<>% mutate(Dataset = if_else(study_id == "Zeng", "Cells", "Nuclei"))

pdf(paste0(base_folder, "cluster_size_binned.pdf"), height=5, width = 7)
for_plot %<>% mutate(n_binned = cut(n, breaks = c(8,20,50,100,1000,10000,1000000)))
for_plot %>% filter(n_binned == "(1e+04,1e+06]")
for_plot %<>% filter(!is.na(n_binned))
levels(for_plot$n_binned)
for_plot %<>% mutate(n_binned = factor(n_binned, labels = c("9-20", "21-50","51-100","101-1,000", "1,001-10,000", "10,001-1,000,000")))
for_plot %>% arrange(n)
ggplot(data = for_plot, aes(x = n_binned, fill = Dataset))  +  
  geom_bar(position = "dodge2") + theme_bw() +xlab("Cluster size") + ylab("Cluster count") +
  theme(legend.position = "bottom") + scale_fill_manual(values = c("#4450A2", "#E21E25"))
dev.off()  


top_hits %<>% inner_join(cell_type_sizes %>% select(`Study_ID|Celltype_1` = cell.type, cell_count_celltype_1 = n))
top_hits %<>% inner_join(cell_type_sizes %>% select(`Study_ID|Celltype_2` = cell.type, cell_count_celltype_2 = n))

correlation <- cor(top_hits  %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(cell_count_celltype_1), top_hits  %>% filter(Match_type == "Reciprocal_top_hit") %>% pull(cell_count_celltype_2), m='s')

#fix scales 
pdf(paste0(base_folder, "recip_cluster_size_correlation.", top_hits %>% nrow(), "_pairs.pdf"), height=6, width = 6)
convert_ticks <- function(x) {
  10**x
}

#scatter plot of recip sizes
ggplot(top_hits %>% filter(Match_type == "Reciprocal_top_hit"), aes(y=log10(cell_count_celltype_2), x = log10(cell_count_celltype_1))) + 
  #geom_point() + theme_bw() + ylab("Nucleus cluster size (log10 scale)") + xlab("Cell cluster size (log10 scale)") + geom_smooth(method = 'lm', se = FALSE) +
  geom_point() + theme_bw() + ylab("Nucleus cluster size") + xlab("Cell cluster size") + geom_smooth(method = 'lm', se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  annotate("text",
           x = 1, 
           y = 5.5, 
           #label = paste0("Spearman r: ", correlation %>% format(digits = 2)),
           label = paste0("r = ", correlation %>% format(digits = 2)),
           hjust = 0, 
           vjust = 1
  ) + scale_x_continuous(labels = convert_ticks) + scale_y_continuous(labels = convert_ticks)
dev.off()
