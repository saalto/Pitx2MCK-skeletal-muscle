## ---- Volcano plot with Reactome Database-identified proteins ----
#' script for making Volcano plots for RNAseq data
#' @param gene_frame: data frame with fold change/p value info
#' @param fc_name: name of column corresponding to log fold change
#' @param p_name: name of column with p-values (adjusted)
#' @param p_levels: a sorted (smallest to largest) vector of three (only 3!) cutoffs for color coding
#' @param group_colors: a set of four colors (1 per p_level +1 for nonsignificant genes)
#' @param file_name: name of the file to write the result
#' @export

volcano_plot2 <- function(gene_frame, gene_frame3,fc_name="logFC",p_name="FDR",p_levels=c(0.05),group_colors=c("#009E73","#000000"),file_name){
  # "#CC79A7","#999999"
  require(ggplot2)
  # append a column to gene_frame that computes -log(p)
  gene_frame[["neglogp"]] <- -1.0*log(gene_frame[[p_name]])
  # create a color column with three levels based on FDR
  # d marks nonsignificant genes
  # c marks genes >= level 1
  # b marks genes between levels 1 and 2
  # a marks genes between levels 2 and 3
  c_vec <- rep("d",dim(gene_frame)[1])
  # c_vec[gene_frame[[p_name]] <= p_levels[3]] <- "c"
  # c_vec[gene_frame[[p_name]] <= p_levels[2]] <- "b"
  c_vec[gene_frame[[p_name]] <= p_levels[1]] <- "a"
  # bind to the data frame
  gene_frame[["cutoff"]] <- c_vec
  # get overall x limits
  max_x <- max(abs(gene_frame[[fc_name]]))
  # get the biggest neglogp
  max_y <- max(gene_frame[["neglogp"]])
  
  # subset of gene frame where genes are present in Reactome database search
  #gene_frame$OutlierNames7 = row.names(gene_frame)
  #idx7 <- as.character(group)
  #gene_frame$OutlierNames7[!(gene_frame$OutlierNames7 %in% idx7)] = NA
  
  # Genes with high DE in TA and low in So (Blue for down-regulated)
  # subset of gene frame with genes with log(fold change) less than -1
  genesDOWN <- gene_frame3[gene_frame3$logFC<(-1),]
  gene_frame$OutlierNames1 = row.names(gene_frame)
  idx1 <- as.character(row.names(genesDOWN))
  gene_frame$OutlierNames1[!(gene_frame$OutlierNames1 %in% idx1)] = NA
  
  # Genes with high DE in So and low in TA (Red for up-regulated)
  # subset of gene frame with genes with a log(fold change) greater than 1
  genesUP <- gene_frame3[gene_frame3$logFC>(1),]
  gene_frame$OutlierNames3 = row.names(gene_frame)
  idx3 <- as.character(row.names(genesUP))
  gene_frame$OutlierNames3[!(gene_frame$OutlierNames3 %in% idx3)] = NA
  
  # create the figure
  pdf(file_name)
  pl <- ggplot2::ggplot(gene_frame, ggplot2::aes_string(x=fc_name,y="neglogp",colour="cutoff")) +
    ggplot2::geom_point(size=1.0,alpha=0.1) + 
    ggplot2::xlim(-max_x,max_x) + 
    ggplot2::ylim(0,1.01*max_y) +
    # fix the colors
    ggplot2::scale_colour_manual(labels=c(paste("p <",p_levels[1]),paste("p >",p_levels[1])),values=group_colors)+
    # remove gray background
    ggplot2::theme_bw()+
    # horizontal line at FDR = 0.05
    # geom_hline(yintercept=-1.0*log(0.05),color="black",linetype="dashed",size=1)+
    
    # font sizes
    # axis, ticks, and legend labels
    ggplot2::theme(axis.title.x = ggplot2::element_text(size=18),
                   axis.title.y = ggplot2::element_text(size=18),
                   axis.text.x = ggplot2::element_text(size=14),
                   axis.text.y = ggplot2::element_text(size=14),
                   legend.title = ggplot2::element_text(size=16),
                   legend.text = ggplot2::element_text(size=12),
                   legend.position="top") +
    # change x and y axis label
    ggplot2::labs(x = expression(paste('log'['2']*'(Fold Change)')), 
                  y = expression(paste('-ln(FDR)'))) +
    
    # override alpha and marker size for legend
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1,size=6))) +
    
    # label genes found in Reactome database search
    # ggplot2::geom_text(ggplot2::aes(label=OutlierNames7), size=1, colour="orange", vjust=0, hjust=-0.1)+
    # label genes found in Reactome database search
    ggplot2::geom_text(ggplot2::aes(label=OutlierNames1), size=1, colour="blue", vjust=0, hjust=-0.1)+  
    # label genes found in Reactome database search
    ggplot2::geom_text(ggplot2::aes(label=OutlierNames3), size=1, colour="red", vjust=0, hjust=-0.1)
    
  print(pl)
  dev.off()
}

## ---- end ----