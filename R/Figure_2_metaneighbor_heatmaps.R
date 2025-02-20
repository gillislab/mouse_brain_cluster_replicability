library(circlize)
library(plotROC)
library(viridis)
library(zoo)
library(dendextend)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(readxl)
library(here)
library(magrittr)
library(MetaNeighbor)
library(reshape2)
library(gplots)
library(conflicted)
library(ComplexHeatmap)
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
all_AUROCs_matrix <- all_AUROCs_matrix[sort(rownames(all_AUROCs_matrix)), sort(rownames(all_AUROCs_matrix))]
dim(all_AUROCs_matrix)

study_colors <- tibble(study = getStudyId(rownames(all_AUROCs_matrix)))

study_colors %<>% mutate(color = recode(study, Macosko = '#E21E25', Zeng = '#4450A2'))
study_colors %<>% mutate(study_type = recode(study, Macosko = "Nuclei", Zeng = "Cells"))
study_colors %>% group_by(study) %>% count()

division_colors <- tibble(study = getStudyId(rownames(all_AUROCs_matrix)), cell_type = getCellType(rownames(all_AUROCs_matrix)))
division_colors %<>% mutate(prefix = gsub("_.*","",cell_type))
division_colors %<>% mutate(prefix = if_else(study == "Zeng", "Zeng", prefix))
division_colors %>% group_by(prefix) %>% count()
division_colors %<>% mutate(prefix = as.factor(prefix))
#default ggplot2 colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
division_colors %<>% mutate(color = as.character(factor(prefix, levels = levels(prefix), labels = gg_color_hue(length(levels(prefix))))))



dendro <- stats::as.dendrogram(orderCellTypes(all_AUROCs_matrix))
dendro_order <- order.dendrogram(dendro)
dendro_table <- tibble(cell_type = rownames(all_AUROCs_matrix)[dendro_order])
dendro_table %<>% mutate(dendro_order = row_number())
dendro_table %<>% mutate(cell_type_short = getCellType(cell_type))
dendro_table %<>% mutate(study = getStudyId(cell_type))
study_colors %<>% mutate(color = recode(study, Macosko = '#E21E25', Zeng = '#4450A2'))
dendro_table %<>% mutate(study_type = recode(study, Macosko = "Nuclei", Zeng = "Cells", Siletti = "Nuclei"))
dendro_table %>% write_csv(paste0(base_folder, "aurocs_dendro_order.csv"))

focus_order <- dendro_table %>% filter(dendro_order < 247)
focus_order %<>% arrange(dendro_order)
focus_order %>% pull(cell_type)
focus_AUROC_matrix <- all_AUROCs_matrix[focus_order %>% arrange(-dendro_order) %>% pull(cell_type), focus_order %>% arrange(dendro_order) %>% pull(cell_type)]
focus_AUROC_matrix[1:5,1:5]

colors = c(seq(0,0.79, length.out=80),seq(0.8,0.89, length.out=10),seq(0.9,1, length.out=11))
my_palette <- c(colorRampPalette(c("#313695", "#74ADD1", "#E0F3F8"))(n = 80),
                colorRampPalette(c("#FFFFBF", "#FEE090", "#FDAE61"))(n = 10),
                colorRampPalette(c("#F46D43", "#D73027", "#A50026"))(n = 11))
color_mapping <- colorRamp2(colors, my_palette)
color_vector <- rep(colors, each = 10)
barplot(rep(1, length(color_vector)), col = color_mapping(color_vector), border = NA, axes = FALSE)

color_mapping <- colorRampPalette(my_palette)

subclass_info <- read_csv(here("data","whole_mouse_brain", "zeng", "from_API", "Cluster_colors_and_tree.abc_atlas_access.csv"))
subclass_info %<>% mutate(class = gsub("\\d+ ", "", class))
subclass_info %<>% mutate(cell_type = paste0("Zeng|", cluster_alias))
subclass_info %<>% select(cell_type, class, class_color)

focus_order %<>% left_join(subclass_info)

subclass_info %>% select(class, class_color) %>% distinct() %>% unite("combined", class, class_color, sep = " ") %>% setNames(., .) -> result

class_col_vec <- subclass_info %>% select(class, class_color) %>% distinct() %>% pull(class_color)
names(class_col_vec )<- subclass_info %>% select(class, class_color) %>% distinct() %>% pull(class)
class_col_vec <- c(class_col_vec, Nuclei = "#FFFFFF")
focus_order %<>% mutate(class = if_else(is.na(class), "Nuclei", class))
focus_order %<>% mutate(class_color = if_else(class == "Nuclei", "#FFFFFF", class_color))

bottom_annotation <- HeatmapAnnotation(df=focus_order %>% arrange(-dendro_order) %>% select(Dataset = study_type) %>% as.data.frame(), 
                                       col = list(Dataset = c("Nuclei"="#E21E25", "Cells"="#4450A2") ), show_annotation_name=FALSE)
left_annotation <- rowAnnotation(df=focus_order %>% arrange(-dendro_order) %>% select(Class = class, Dataset = study_type) %>% as.data.frame(),
                                 col = list(Dataset = c("Nuclei"="#E21E25", "Cells"="#4450A2"), 
                                            Class = class_col_vec), show_annotation_name=FALSE)#, na_col="white")

pdf(paste0(base_folder, "zoomed_auroc_corner_scale_plus_class.pdf"), height=6, width = 8)
focus_AUROC_matrix_rev <- focus_AUROC_matrix[, rev(colnames(focus_AUROC_matrix))]
Heatmap(focus_AUROC_matrix_rev, 
        show_row_dend= FALSE, 
        show_column_dend = FALSE, 
        show_row_names = FALSE, 
        show_column_names = FALSE,
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        col = color_mapping, 
        bottom_annotation = bottom_annotation, 
        left_annotation = left_annotation, 
        name="AUROC")
dev.off()


#################
#write out color bar for the whole dendrogram/heatmap
#################

png(paste0(base_folder, "aurocs_all.png"), width=14000, height = 14000)
plotHeatmap(all_AUROCs_matrix, cex = 0.5,  labRow = FALSE, labCol = FALSE, labCol = FALSE,
            ColSideColors = study_colors %>% pull(color),
            RowSideColors = study_colors %>% pull(color), 
            col = color_mapping)
dev.off() #this was flipped on X for the figure


dendro_table %<>% left_join(subclass_info)
dendro_table
color_table <- dendro_table %>% select(class, class_color) %>% distinct()

pdf(paste0(base_folder, "aurocs_all_heatmap_class_color_bar.pdf"), height=30, width = 4)
ggplot(dendro_table, aes(x = 1, y = dendro_order, fill = class)) +
  geom_tile() +
  scale_fill_manual("Color", values = setNames(color_table %>% pull(class_color), color_table %>% pull(class)), na.value="white") +
  theme_minimal() + ylab("") + xlab("") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())
dev.off()

#write out the study color bar, final figure is assembled in Illustrator
study_colors %<>% mutate(color = recode(study, Macosko = '#E21E25', Zeng = '#4450A2'))
dendro_table %<>% mutate(study_type = recode(study, Macosko = "Nuclei", Zeng = "Cells"))
dendro_table %<>% mutate(study_color = if_else(study == "Zeng",'#4450A2', '#E21E25'))
dendro_table  %>% select(study_type, study_color) %>% distinct()
pdf(paste0(base_folder, "aurocs_all_heatmap_dataset_color_bar.pdf"), height=30, width = 4)
ggplot(dendro_table, aes(x = 1, y = dendro_order, fill = study_type)) +
  geom_tile() +
  scale_fill_manual("Dataset", values = setNames(dendro_table %>% pull(study_color), dendro_table %>% pull(study_type)), na.value="white") +
  theme_minimal() + ylab("") + xlab("") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())
dev.off()



