#' script for assembling GO data for use with GOplot. returns a data frame of the circ
#'  type (see GOplot docs) and a matrix for use with GOplot (see chord_dat docs)
#' @param gene_frame: data frame with fold change/p value info, as well as gene names
#' @param goplot_data: list containing results from GO_assemble script
#' @param file_name: file name of the resulting chord plot
#' @export

scr_chordplot2 <- function(gene_frame, goplot_data, file_name){
  require(GOplot)

  # ----- SHOULD NOT NEED TO CHANGE STUFF DOWN HERE ----
  # set overall color scale (but make it symmetric) using max/min in table
  gene_frame <- data.frame(gene_frame)
  max_fc <- max(gene_frame$logFC)
  min_fc <- min(gene_frame$logFC)
  if (max_fc > abs(min_fc)){
    fc_max <- max_fc
    fc_min <- -1.0*max_fc
  } else {
    fc_max <- -1.0*min_fc
    fc_min <- min_fc
  }

  # make the plot
  pdf(file_name, width = 20, height = 20)
  pl <- GOChord(goplot_data, space=0.02,gene.order='logFC',gene.space=0.25,lfc.col=c('red','white','blue'),lfc.min=fc_min,lfc.max=fc_max,border.size=0, process.label = 10)
  pl <- pl + coord_fixed()
  print(pl)
  dev.off()
}
