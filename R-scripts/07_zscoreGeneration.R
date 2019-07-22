# ---- Section 1: Heat map of Soleus data ----
  # - Creating the Read Counts variable -
  AllSoleus_ReadCounts <- read.delim(file="./data/soleus/soleus-data-combined.txt", header = FALSE, sep=",")
  AllSoleus_ReadCounts1 <- cbind(C1= AllSoleus_ReadCounts$V2, C2= AllSoleus_ReadCounts$V3, C3= AllSoleus_ReadCounts$V4, 
                                 C4= AllSoleus_ReadCounts$V5, C5= AllSoleus_ReadCounts$V6, C6= AllSoleus_ReadCounts$V7, 
                                 M1= AllSoleus_ReadCounts$V8, M2= AllSoleus_ReadCounts$V9, M3= AllSoleus_ReadCounts$V10, 
                                 M4= AllSoleus_ReadCounts$V11, M5= AllSoleus_ReadCounts$V12)
  row.names(AllSoleus_ReadCounts1) <- AllSoleus_ReadCounts$V1
  
  # - Calculating the z-scores -
  AllSoleus_DEMetGenes_ReadCounts <- AllSoleus_ReadCounts1[row.names(AllSoleus_ReadCounts1) %in% deedgSOL_flt$sig_genes,]
  df1 <- AllSoleus_DEMetGenes_ReadCounts[apply(AllSoleus_DEMetGenes_ReadCounts,1,prod)>0,]
  df2 <- log(df1)
  AllSoleus_zscore <- apply(df2,1,scale)
  
  # - Creating the heatmap -
  breakList = seq(min(AllSoleus_zscore),max(AllSoleus_zscore), by=0.2)
  colorz = colorRampPalette(colorz)(length(breakList))
  
  pheatmap(AllSoleus_zscore, breaks=breakList, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
           show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0, filename = "SoMutantvsSoControl_heatmap.jpeg")

# ---- end ----
  
# ---- Section 2: Heat map of Soleus data set without Sample C6 ----
  # Sample C6 was removed due to low z-scores
  # - Creating the Read Counts variable -
  AllSoleusWO6_ReadCounts <- read.delim(file="./data/soleus/modified/soleus-WO6-ReadCounts.txt", header = FALSE, sep = ",")
  AllSoleusWO6_ReadCounts1 <- cbind(C1= AllSoleusWO6_ReadCounts$V2, C2= AllSoleusWO6_ReadCounts$V3, 
                                    C3= AllSoleusWO6_ReadCounts$V4, C4= AllSoleusWO6_ReadCounts$V5, 
                                    C5= AllSoleusWO6_ReadCounts$V6, M1= AllSoleusWO6_ReadCounts$V7, 
                                    M2= AllSoleusWO6_ReadCounts$V8, M3= AllSoleusWO6_ReadCounts$V9, 
                                    M4= AllSoleusWO6_ReadCounts$V10, M5= AllSoleusWO6_ReadCounts$V11)
  row.names(AllSoleusWO6_ReadCounts1) <- AllSoleusWO6_ReadCounts$V1
  
  # - Calculating the z-scores -
  AllSoleusWO6_DEMetGenes_ReadCounts <- AllSoleusWO6_ReadCounts1[row.names(AllSoleusWO6_ReadCounts1) %in% deedgSOLWO6_flt$sig_genes,]
  df1 <- AllSoleusWO6_DEMetGenes_ReadCounts[apply(AllSoleusWO6_DEMetGenes_ReadCounts,1,prod)>0,]
  df2 <- log(df1)
  AllSoleusWO6_zscore <- apply(df2,1,scale)
  
  # - Creating the heatmap -
  # Note: this soleus data set was used to scale the color legend
  # because this data set had the highest and lowest z-scores
  breakList = seq(min(AllSoleusWO6_zscore),max(AllSoleusWO6_zscore), by=0.2)
  colorz = colorRampPalette(colorz)(length(breakList))
  
  pheatmap(AllSoleusWO6_zscore, breaks=breakList, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
           show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0, filename = "SoMutantvsSoControl_WO6_heatmap.jpeg")
# ---- end ----
  
# ---- Section 3: Heat map of Tibialis Anterior data ----
  # - Creating the Read Counts variable -
  AllTibialis_ReadCounts <- read.delim(file = "./data/tibialis anterior/tibialis-data-combined.txt", header=FALSE, sep = ",")
  AllTibialis_ReadCounts1 <- cbind(C1=AllTibialis_ReadCounts$V2, C2=AllTibialis_ReadCounts$V3, C3=AllTibialis_ReadCounts$V4, 
                                   M2=AllTibialis_ReadCounts$V5, M3=AllTibialis_ReadCounts$V6, M4=AllTibialis_ReadCounts$V7, 
                                   M5=AllTibialis_ReadCounts$V8)
  row.names(AllTibialis_ReadCounts1) <- AllTibialis_ReadCounts$V1
  
  # - Calculating the z-scores -
  AllTibialis_DEMet_ReadCounts <- AllTibialis_ReadCounts1[row.names(AllTibialis_ReadCounts1) %in% deedgTA_flt$sig_genes,]
  df1 <- AllTibialis_DEMet_ReadCounts[apply(AllTibialis_DEMet_ReadCounts,1,prod)>0,]
  df2 <- log(df1)
  AllTibialis_zscores <- apply(df2,1,scale)
  
  # calculating the z-scores of soleus data set because that data set had the highest and lowest z-scores
  # those values were used to scale the colors uniformly across both soleus and tibialis data sets
  # the data set using was without Sample C6 due to low z-scores
  
  # Creating the Read Counts variable -
  AllSoleusWO6_ReadCounts <- read.delim(file="./data/soleus/modified/soleus-WO6-ReadCounts.txt", header = FALSE, sep = ",")
  AllSoleusWO6_ReadCounts1 <- cbind(C1= AllSoleusWO6_ReadCounts$V2, C2= AllSoleusWO6_ReadCounts$V3, 
                                    C3= AllSoleusWO6_ReadCounts$V4, C4= AllSoleusWO6_ReadCounts$V5, 
                                    C5= AllSoleusWO6_ReadCounts$V6, M1= AllSoleusWO6_ReadCounts$V7, 
                                    M2= AllSoleusWO6_ReadCounts$V8, M3= AllSoleusWO6_ReadCounts$V9, 
                                    M4= AllSoleusWO6_ReadCounts$V10, M5= AllSoleusWO6_ReadCounts$V11)
  row.names(AllSoleusWO6_ReadCounts1) <- AllSoleusWO6_ReadCounts$V1
  
  # Calculating the z-scores
  AllSoleusWO6_DEMetGenes_ReadCounts <- AllSoleusWO6_ReadCounts1[row.names(AllSoleusWO6_ReadCounts1) %in% deedgSOLWO6_flt$sig_genes,]
  df1 <- AllSoleusWO6_DEMetGenes_ReadCounts[apply(AllSoleusWO6_DEMetGenes_ReadCounts,1,prod)>0,]
  df2 <- log(df1)
  AllSoleusWO6_zscore <- apply(df2,1,scale)
  
  # color legend scaling
  breakList = seq(min(AllSoleusWO6_zscore),max(AllSoleusWO6_zscore), by=0.2)
  colorz = colorRampPalette(colorz)(length(breakList))
  
  # - Creating the heatmap -
  pheatmap(AllTibialis_zscores, breaks = breakList,color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
           show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0, filename = "TAMutantvsTAControl_heatmap.jpeg")
# ---- end ----
  
# ---- Section 4: Heat map of TA data set with renamed files ----
  # - Creating the Read Counts variable -
  AllTibialis_renamed_ReadCounts <- read.delim(file = "./data/tibialis anterior/renamed_files/AllTibialis_renamed_ReadCounts.txt", header=FALSE, sep = ",")
  AllTibialis_renamed_ReadCounts1 <- cbind(C1=AllTibialis_renamed_ReadCounts$V2, C2=AllTibialis_renamed_ReadCounts$V3, 
                                  M1=AllTibialis_renamed_ReadCounts$V4, M2=AllTibialis_renamed_ReadCounts$V5, 
                                  M3=AllTibialis_renamed_ReadCounts$V6, M4=AllTibialis_renamed_ReadCounts$V7, 
                                  M5=AllTibialis_renamed_ReadCounts$V8)
  row.names(AllTibialis_renamed_ReadCounts1) <- AllTibialis_renamed_ReadCounts$V1
  
  # - Calculating the z-scores -
  AllTibialis_renamedDEMet_ReadCounts <- AllTibialis_renamed_ReadCounts1[row.names(AllTibialis_renamed_ReadCounts1) %in% deedgTA_flt_renamed$sig_genes,]
  df1 <- AllTibialis_renamedDEMet_ReadCounts[apply(AllTibialis_renamedDEMet_ReadCounts,1,prod)>0,]
  df2 <- log(df1)
  AllTibialis_renamed_zscores <- apply(df2,1,scale)
  
  # - Creating the heatmap -
  # NOTE: The color scale was set based on Soleus data because that data set had the highest and lowest z-scores
  pheatmap(AllTibialis_renamed_zscores, breaks = breakList,color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
           show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0, filename = "TAMutantvsTAControl_renamed_heatmap.jpeg")
# ---- end ----