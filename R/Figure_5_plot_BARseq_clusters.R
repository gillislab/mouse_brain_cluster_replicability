library(ggplot2)
library(SingleCellExperiment)
library(magrittr)
library(dplyr)
library(S4Vectors)
library(here)
library(readr)
library(tibble)
library(conflicted)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)

base_barseq_folder <- here("data/whole_mouse_brain/processed/barseq/")                           

#load labels from https://data.mendeley.com/datasets/8bhhk7c5n9/1
labels <- read_csv(paste0(base_barseq_folder, "labels_20211208_withgaba_nonexc.csv"))
labels %>% pull(cluster) %>% unique() %>% sort() %>% length() #127 after removing NA value

#propigate labels down to cluster level for NA's
labels %<>% mutate(subclass = if_else(is.na(subclass), class, subclass))
labels %<>% mutate(cluster = if_else(is.na(cluster), subclass, cluster))
labels %>% pull(cluster) %>% unique() %>% sort()
labels %>% pull(cluster) %>% unique() %>% sort() %>% length()

#load rds from /data/fischer/ct_assessment/neuron_position.csv
bar_seq_data <- readRDS(paste0(base_barseq_folder, "barseq_210630.rds"))

#join labels
col_data <- colData(bar_seq_data) %>% as.data.frame() %>% rownames_to_column(var = "row_name") %>% tibble()
col_data %<>% mutate(sample = paste0(slice, "_", row_name))
col_data %<>% select(sample, everything())
col_data %<>% left_join(labels %>% select(-notes))
col_data %>% group_by(is.na(cluster)) %>% count()
col_data %>% group_by(is.na(cluster), is.na(ccf_x)) %>% count() #Cross expression also filtered on CCF data

col_data %>% filter(!is.na(cluster), !is.na(ccf_x))  %>% pull(cluster) %>% unique() %>% length()

bar_seq_data <- NULL #to save memory

new_cords <- col_data %>% select(slice, cluster_H3 = cluster, pos_x, pos_y)

#uncomment/comment the below lines to select a cluster - ran interactivley
target <- "P D?"
slice_id <- 2

# target <- "SubCTX_7"
# slice_id <- 23
# 
# target <- "SubCTX_4"
# slice_id <- 7


#will show number of target cells across slices
new_cords %>% filter(cluster_H3 == target) %>% group_by(slice) %>% count() %>% arrange(-n)

#some slices need rotating, this isn't a full check of them all
y_multiplier <- -1
x_multiplier <- 1
if (slice_id %in% c(2, 21, 23, 34, 19, 10,32,4,20,9,3,24)) {
  y_multiplier <- 1
  x_multiplier <- -1
}

new_cords %>% filter(cluster_H3 == target) %>% group_by(slice) %>% count() %>% arrange(slice) %>% as.data.frame()
new_cords %<>% filter(!is.na(cluster_H3))

#uncomment for PDF files
#pdf(paste0(base_folder, "BARseq_", target, "_slice_", slice_id, ".pdf"), height=5, width = 7)
for_plot <- new_cords %>%
  filter(slice == slice_id) %>% mutate(is_target_cluster = cluster_H3 == target)
hits_on_slice <- for_plot %>% group_by(is_target_cluster) %>% dplyr::count() %>% filter(is_target_cluster) %>% pull(n)
point_size <- 0.15
ggplot(data = for_plot, aes(x = pos_y * x_multiplier, y = y_multiplier*pos_x, col = is_target_cluster)) + geom_point(size = 0.15, color ="gray90") +
  geom_point(data = for_plot %>% filter(is_target_cluster), size = point_size, color = "black") + theme_classic() +
  scale_color_discrete(na.value="gray90") + coord_equal() + theme(legend.position="none") +
  #ggtitle(paste0(target, " slice:", slice_id)) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(
    axis.line = element_line(color = 'white'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  )
#dev.off()

#Code to plot all slices - orientations are different up in some slices
for_plot <- new_cords %>% mutate(is_target_cluster = cluster_H3 == target)
hits_on_slice <- for_plot %>% group_by(is_target_cluster,slice) %>% dplyr::count() %>% filter(is_target_cluster) %>% pull(n)
ggplot(data = for_plot, aes(x = pos_y, y = -1*pos_x, col = is_target_cluster)) + 
  geom_point(size = 0.02, color ="gray90") + 
  geom_point(data = for_plot %>% filter(is_target_cluster), size = 0.02, color = "black") + theme_classic() +
  scale_color_discrete(na.value="gray90") + coord_equal() + theme(legend.position="none") +
  facet_wrap(. ~ slice)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(
    axis.line = element_line(color = 'white'),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  ) 
