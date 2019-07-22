#' script for assembling GO data for use with GOplot. returns a data frame of the circ
#'  type (see GOplot docs) and a matrix for use with GOplot (see chord_dat docs)
#' @param gene_frame: data frame with fold change/p value info, as well as gene names
#' @param all_go: data frame with the ENTREZ ID,GO ID, GO evidence, Ontology, and gene name
#' @param double_hash: symbol-to-entrez and entrez-to-symbol hash map
#' @param overrep_file: Over-represented GO term list generated from catenrich.R function
#' @param go_ids: vector of IDs to check against
#' @param annot: organism database
#' @export
GO_assemble <- function(gene_frame,double_hash,all_go, overrep_file,go_ids, annot=org.Mm.eg.db){
  require(GOplot)
  # these will hold the results as we loop through IDs
  ID = c()
  term = c()
  count = c()
  genes = c()
  logFC = c()
  
  # go through GO ids one-by-one and pull gene information
  for (i in 1:length(go_ids)){
    go_id2 = go_ids[i]
    # pull the gene list for specific GO ID
    g <- all_go[all_go$GOALL %in% go_id2,]$GENES
    genes <- c(genes,g)
    # now info
    count <- c(count,rep(length(genes),length(g)))
    ID <- c(ID,rep(go_id2,length(g)))
  }
  
  # category is BP when dealing with GO ids
  category <- rep('BP',length(genes))
  
  # pull logFC and term
  for(i in 1:length(genes)){
    g = genes[i]
    ident = ID[i]
    logFC <- c(logFC, gene_frame[row.names(gene_frame) %in% g,])
    overrep_term <- overrep_file$overrep[overrep_file$overrep$category %in% ident,]
    term <- c(term, overrep_term$term)
  }
  
  # this mimics the data frame produced by circ_dat
  circ <- data.frame(category,ID,count,term,genes,logFC,stringsAsFactors=FALSE)
  # produce the one that chord data make; need the "genes" data frame
  u_genes <- unique(circ$genes)
  u_logFC <- unique(circ$logFC)
  genes_for_chord = data.frame("ID"=u_genes,"logFC"=u_logFC,stringsAsFactors=FALSE)
  #genes_for_chord = genes_for_chord[genes_for_chord$ID %in% genelist$GENE,]
  # this is the one that chord_dat makes
  chord <- chord_dat(circ,genes_for_chord,unique(circ$term))
  # put them both in a list
  goplot_data <- list("circ" = circ, "chord" = chord)
}
