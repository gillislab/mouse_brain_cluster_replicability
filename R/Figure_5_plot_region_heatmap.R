library(corrplot)
library(RColorBrewer)


# regions of interest
regs = c('CB', 'HPF', 'STRd', 'OLF', 'Isocortex', 'BS', 'MB', 'TH', 'PAL', 'HY')
my_palette <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"))(100))

# function to get change in prop for heatmap
get_change_in_prop <- function(df, regs){
    mat1 = df[match(regs, df$region), c('prop_recip_hit_Zeng', 'prop_recip_hit_Macosko')]
    mat2 = colMeans(mat1, na.rm = T)
    
    mat1[,1] <- mat1[,1] - mat2[1]   # get change in prop.
    mat1[,2] <- mat1[,2] - mat2[2]   # get change in prop.
    mat1 <- t(as.matrix(mat1))
    colnames(mat1) = regs
    return(mat1)
}

# load file with number of cells with Zeng/Macosko reciprocal hits per region
# STR region called as STRd in this file
mtd = read.delim('~/results/other_mouse_studies/reciprocal_hits_by_region.csv', sep = ',')


# get change in proportion of no. of reciprocal hit cells with Zeng/Macosko
mat = get_change_in_prop(mtd, regs)
rownames(mat) = c('Zeng', 'Macosko')


# plot heat map
pdf('prop-cell-mapped-region-heatmap.pdf', width = 8, height = 5)
corrplot(mat, method = 'color', is.corr = F, col = my_palette, tl.col = 'black',
        na.label = 'square', na.label.col = '#f6f6f6', col.lim = c(-0.4, 0.4))
dev.off()
