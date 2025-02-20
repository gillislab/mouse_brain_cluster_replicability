library(grDevices)
library(here)
library(ggplot2)
library(reshape2)
library(magrittr)
library(dplyr)
library(MetaMarkers)
library(readr)
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::first)
conflicts_prefer(dplyr::count)

Sys.setenv("VROOM_CONNECTION_SIZE"= 21000000)

greedy_tests_Mac <- read_csv(here("results", "MetaMarkers_2009_reciprocals", "MetaMarkers_greedy_auc_tests.Macosko.1706156897.csv.gz"))
greedy_tests_Zeng <- read_csv(here("results", "MetaMarkers_2009_reciprocals", "MetaMarkers_greedy_auc_tests.Zeng.1706190036.csv.gz"))
threshold_value <- 0.95
thresholded <- greedy_tests_Zeng %>% filter(auroc > threshold_value) %>% 
  group_by(genes_used) %>% count()%>% filter(genes_used < 51) %>% mutate(Study = "Zeng")
thresholded %<>% bind_rows(greedy_tests_Mac %>% filter(auroc > threshold_value) %>% group_by(genes_used) %>% count()%>% filter(genes_used < 51) %>% mutate(Study = "Macosko"))
thresholded %<>% mutate(Dataset = if_else(Study=="Macosko", "Nuclei", "Cells"))
thresholded %>% group_by(Dataset) %>% filter(genes_used==4)
thresholded %>% group_by(Dataset) %>% slice_max(n)


thresholded %<>% mutate(color = recode(Study, Macosko = '#E21E25', Zeng = '#4450A2'))
pdf(here("results", "MetaMarkers_2009_reciprocals", "thresholded_greedy_metamarkers.pdf"), height=5, width = 6)
ggplot(data = thresholded  %>% filter(genes_used < 51), aes(x = genes_used, y = n, color=Dataset)) + geom_line() + theme_bw() + 
  ylab(paste0("Clusters with best-vs-next AUROC > ", threshold_value)) + xlab("Markers used") + scale_color_manual(values = c("Nuclei" = '#E21E25', "Cells" = '#4450A2')) +
  scale_x_continuous(expand = c(0, 0))
dev.off()  
