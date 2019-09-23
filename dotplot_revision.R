library(ggplot2)
library(dplyr)
library(stringr)
## for reordering the factor
library(forcats) 
setSize = 1062

## count the gene number
gene_count<- yy3@result %>% group_by(ID) %>% summarise(count = sum(str_count(geneID, "/")) + 1)

## merge with the original dataframe
dot_df<- left_join(yy3@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

dot_df = dot_df[1:50,] ## small dataset
dot_df$type = "upregulated"
# I do not know where the NES is though
dot_df$type[dot_df$NES < 0] = "downregulated"

## plot
p <- ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")

p + facet_grid(.~type)