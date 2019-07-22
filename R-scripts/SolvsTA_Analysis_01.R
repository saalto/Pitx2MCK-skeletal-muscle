# --- Load the libraries of pheatmap and RColorBrewer---
  library(pheatmap)
  library(RColorBrewer)
  colorz = rev(brewer.pal(10,"RdBu"))
  
  # set variables and source of necessary functions
  GeneUniverse <- readRDS(file="./data/GeneUniverse.RDS")
  mouse_hash <- readRDS(file = "./data/mouse-hash.rds")
  enrich_type <- "GO:BP"
  source("./R-scripts/03_catenrich.R")
  source("./R-scripts/04_GO_genes.R")
  source("./R-scripts/05_GO_assemble.R")
  source("./R-scripts/06_scr_chordplot2.R")
  
# ---- (7) Filter out genes associated with Metabolic Process (GO:0008152) via PANTHER ----
  # filter out genes in soleus data set 
  soleus_WO6_DEGenes <- deedgSOLWO6_flt$table[row.names(deedgSOLWO6_flt$table) 
                                              %in% deedgSOLWO6_flt$sig_genes,]
  soleus_WO6_DEMetGenes <- soleus_WO6_DEGenes[row.names(soleus_WO6_DEGenes) %in% GeneUniverse, ]
  # filter out genes in tibialis data set
  tibialis_renamed_DEGenes <- deedgTA_flt_renamed$table[row.names(deedgTA_flt_renamed$table) %in% 
                                                          deedgTA_flt_renamed$sig_genes,]
  tibialis_renamed_DEMetGenes <- tibialis_renamed_DEGenes[row.names(tibialis_renamed_DEGenes) 
                                                          %in% GeneUniverse,]
# ---- end ----
  
# ---- (8) analyze GO term enrichment of soleus and tibialis separately----
  # run the following
  soleus_WO6_enrich <- catenrich(all.genes = GeneUniverse, sig.genes = row.names(soleus_WO6_DEMetGenes), 
                                 entrez_hash = mouse_hash, enrich_type = enrich_type)
  tibialis_renamed_enrich <- catenrich(all.genes = GeneUniverse, sig.genes = row.names(tibialis_renamed_DEMetGenes), 
                                       entrez_hash = mouse_hash, enrich_type = enrich_type)
  # identify shared GO terms
  SharedSOLvsTA <- Reduce(intersect, list(soleus_WO6_enrich$overrep$category, 
                                          tibialis_renamed_enrich$overrep$category))
  # load the genes associated with those shared GO terms
  soleus_sig_genes <- readRDS(file="./data/soleus_overrep_genes.RDS")
  tibialis_sig_genes <- readRDS(file="./data/tibialis_overrep_genes.RDS")
  # concatenate the list of genes
  soleus_sig_genes2 <- Reduce(union, soleus_sig_genes)
  tibialis_sig_genes2 <- Reduce(union, tibialis_sig_genes)
  # identify shared genes between the two data sets
  SharedSOLvsTA2 <- Reduce(intersect, list(soleus_sig_genes2, tibialis_sig_genes2))
# ---- end ----
  
# ---- (9) calculate the average log-fold expression combining soleus and tibialis anterior data sets ----
  # create a table of the log-fold changes
  # pull information from soleus data set
  soleus_joined_sig_genes <- deedgSOLWO6_flt$table[row.names(deedgSOLWO6_flt$table) %in% SharedSOLvsTA2,]
  soleus_joined_sig_genes1 <- soleus_joined_sig_genes[order(row.names(soleus_joined_sig_genes)),]
  # pull information from tibialis data set
  tibialis_joined_sig_genes <- deedgTA_flt_renamed$table[row.names(deedgTA_flt_renamed$table) %in% SharedSOLvsTA2,]
  tibialis_joined_sig_genes1 <- tibialis_joined_sig_genes[order(row.names(tibialis_joined_sig_genes)),]
  
  # combine and make new data frame
  sharedFC <- cbind(Soleus_FC= soleus_joined_sig_genes1$logFC, Tibialis_FC= tibialis_joined_sig_genes1$logFC)
  row.names(sharedFC) <- row.names(soleus_joined_sig_genes1)
  sharedFC <- data.frame(sharedFC)
  # take the average of genes with same log-fold change sign between both data sets
  averages <- rowMeans(sharedFC)
  sharedFC <- cbind(sharedFC, Average=averages)
  
  # take the absolute value of the average to remove log-fold change less than 1
  removal <- abs(sharedFC$Average)
  sharedFC <- cbind(sharedFC, Absolute=removal)
  # compare log-fold change signs
  checkSign <- sign(sharedFC$Soleus_FC)*sign(sharedFC$Tibialis_FC)
  sharedFC <- cbind(sharedFC, Sign=checkSign)
  sharedFC <- sharedFC[with(sharedFC, order(sharedFC$Sign, sharedFC$Average, decreasing = TRUE)),]
  # remove opposing log-fold change signs
  sharedFC1 <- sharedFC[sharedFC$Sign == 1,]
  # remove average log-fold change less than the absolute value of one
  sharedFC1 <- sharedFC1[!sharedFC1$Absolute < 1,]
# ---- end ----
  
# ---- (10) generate heatmap of the shared genes and their log-fold changes ----
  # isolate the greater than 1 log-fold change of shared metabolic process between soleus and tibialis 
  sharedFC2 <- cbind(Soleus=sharedFC1$Soleus_FC, Tibialis=sharedFC1$Tibialis_FC)
  row.names(sharedFC2) <- row.names(sharedFC1)
  # creating Figure 3 heatmap of the log-fold changes
  breakList = seq(min(sharedFC2),max(sharedFC2), by=0.2)
  colorz = colorRampPalette(colorz)(length(breakList))
  pheatmap(sharedFC2, breaks=breakList, color = colorz, cluster_rows = FALSE, show_colnames = TRUE, 
           show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0, filename = "./SovsTA_DEMet_heatmap.jpeg")
# ---- end ----
  
# ---- (11) identify GO terms based on ENTREZ IDs and GO evidence ----
  # run the following:
  sharedSOLvsTA3 <- cbind(logFC=sharedFC1$Average)
  row.names(sharedSOLvsTA3) <- row.names(sharedFC1)
  sharedSOLvsTA3 <- data.frame(sharedSOLvsTA3)
  sharedSOLvsTA_allgo <- GO_genes(gene_frame = sharedSOLvsTA3, double_hash = mouse_hash, 
                                  go_ids = SharedSOLvsTA)
  sharedSOLvsTA_allgo <- sharedSOLvsTA_allgo[sharedSOLvsTA_allgo$GENES %in% row.names(sharedSOLvsTA3),]
# ---- end ----
  
# ---- (12) create variable to make chord plots of GO terms and corresponding genes ----
  # run the following:
  sharedSOLvsTA_goplot <- GO_assemble(gene_frame = sharedSOLvsTA3, double_hash = mouse_hash,
                                    all_go = sharedSOLvsTA_allgo, overrep_file = soleus_WO6_enrich,
                                    go_ids = SharedSOLvsTA)
  # set universal variables
  group_total <- sharedSOLvsTA_goplot$chord[ , colSums(sharedSOLvsTA_goplot$chord !=0)>0]
  group_logFC <- c("logFC")
  group_logFC2 <- group_total[ , (colnames(group_total) %in% group_logFC)]
# ---- end ----
  
# ---- (13) generate chord plots of specific GO terms----
  file.edit(title = "./R-scripts/SolvsTA_Analysis_02.R")

# ---- end ----
  