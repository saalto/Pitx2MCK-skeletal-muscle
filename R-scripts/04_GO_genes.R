#' Returns a list of genes, taken from an input list, that have a given GO term ID
#' @param gene_frame: list of gene symbols
#' @param double_hash: R hash() giving symbol-to-entrez ($stoe) mapping AND entrez-to-symbol ($etos) mapping
#' @param go_ids: go_id to search for
#' @param annot: AnnotationDBI object to use
#' @param specific: FALSE results in generic terms and is longer to process
#' @export

GO_genes <- function(gene_frame,double_hash,go_ids,annot=org.Mm.eg.db, specific=FALSE){
  require(hash)
  require(org.Mm.eg.db)
  # because of missing ENTREZ IDs, some requested genes in gene_list may not be matchable.
  # we therefore create a modified gene_list without these
  mod_genes <- gene_frame[unname(has.key(row.names(gene_frame),double_hash$stoe)),]
  # use the hash to get entrez ids for the gene list
  entrez_ids <- hash::values(double_hash$stoe,keys=row.names(mod_genes),USE.NAMES=FALSE)
  # now get a big data frame with all the GO IDs; get more general terms if specific == FALSE
  if(specific){
    all_go <- AnnotationDbi::select(annot,keys=entrez_ids,columns="GO")
  } else {
    all_go <- AnnotationDbi::select(annot,keys=entrez_ids,columns="GOALL", keytype="ENTREZID")
  }
  # subset the data frame to:
  # - ignore any terms whose evidence comes from ND, NR, or IEA
  all_go <- all_go[all_go$EVIDENCEALL != "ND" & all_go$EVIDENCEALL != "NR" & all_go$EVIDENCEALL != "IEA",]
  # - remove NAs from the result
  all_go <- all_go[!is.na(all_go$GO),]
  # - filter the frame to remove anything other than BP we want
  all_go <- all_go[all_go$ONTOLOGY == "BP",]
  # - filter the frame to remove anything other than the GO term we want
  all_go <- all_go[all_go$GO %in% go_ids,]
  
  # now use the entrez_ids to get back the gene names
  gene_syms <- hash::values(double_hash$etos,keys=all_go$ENTREZID,USE.NAMES=FALSE)
  all_go$GENES <- gene_syms
  return(all_go)
}
