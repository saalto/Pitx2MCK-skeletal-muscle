# ---- Section 1: Heat map of Soleus and Tibialis Anterior data ----
require(pheatmap)
require(RColorBrewer)
colorz = rev(brewer.pal(10,"PiYG"))

# - Creating the Read Counts variable -
SolvsTA_control_ReadCounts <- read.delim(file="c:/Users/sarah/OneDrive/Documents/2018/03_2018_Summer/iteration2/control/TAvsSOL/original files/soleusandTAC-data-combined.txt",
                                         header = FALSE, sep=",")
SolvsTA_control_ReadCounts1 <- cbind(Sol1= SolvsTA_control_ReadCounts$V2, Sol2= SolvsTA_control_ReadCounts$V3, 
                               Sol3= SolvsTA_control_ReadCounts$V4, Sol4= SolvsTA_control_ReadCounts$V5, 
                               Sol5= SolvsTA_control_ReadCounts$V6, TA1= SolvsTA_control_ReadCounts$V7, 
                               TA2= SolvsTA_control_ReadCounts$V8)
row.names(SolvsTA_control_ReadCounts1) <- SolvsTA_control_ReadCounts$V1

# filtering genes from >24421 to subset of statistical significant DE genes (6459 genes)
SolvsTA_control_DEGenes_ReadCounts <- SolvsTA_control_ReadCounts1[row.names(SolvsTA_control_ReadCounts1) 
                                                                     %in% deedg_SOLvsTA_control$sig_genes,]
# Transform read counts to z-scores
df1 <- SolvsTA_control_DEGenes_ReadCounts[apply(SolvsTA_control_DEGenes_ReadCounts,1,prod)>0,] 
df2 <- log(df1) 
SolvsTA_control_zscore <- apply(df2,1,scale)
# 6123 genes present in this list

# label the rows
row.names(SolvsTA_control_zscore) <- c("TA 1", "TA 2", "So 1", "So 2", "So 3", "So 4", "So 5")
# reorder the rows for the image
SolvsTA_control_zscore2 <- data.frame(SolvsTA_control_zscore)
SolvsTA_control_zscore3 <- rbind("So 1"=SolvsTA_control_zscore2["So 1",],
                                 "So 2"=SolvsTA_control_zscore2["So 2",],
                                 "So 3"=SolvsTA_control_zscore2["So 3",],
                                 "So 4"=SolvsTA_control_zscore2["So 4",],
                                 "So 5"=SolvsTA_control_zscore2["So 5",],
                                 "TA 1"=SolvsTA_control_zscore2["TA 1",],
                                 "TA 2"=SolvsTA_control_zscore2["TA 2",]) 
# - Creating the heatmap -
breakList = seq(min(SolvsTA_control_zscore),max(SolvsTA_control_zscore), by=0.2)
#colorz = colorRampPalette(colorz)(length(breakList))
colorz = scico(n=length(breakList), direction = -1, palette = 'roma')
pheatmap(SolvsTA_control_zscore3, breaks=breakList, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0)


# all 24,421 genes
df3 <- SolvsTA_control_ReadCounts1[apply(SolvsTA_control_ReadCounts1,1,prod)>0,]
df4 <- log(df3)
SolvsTA_control_zscore4 <- apply(df4,1,scale)
# label the rows
row.names(SolvsTA_control_zscore4) <- c("TA 1", "TA 2", "So 1", "So 2", "So 3", "So 4", "So 5")
# reorder the rows for the image
SolvsTA_control_zscore5 <- data.frame(SolvsTA_control_zscore4)
SolvsTA_control_zscore6 <- rbind("So 1"=SolvsTA_control_zscore5["So 1",],
                                 "So 2"=SolvsTA_control_zscore5["So 2",],
                                 "So 3"=SolvsTA_control_zscore5["So 3",],
                                 "So 4"=SolvsTA_control_zscore5["So 4",],
                                 "So 5"=SolvsTA_control_zscore5["So 5",],
                                 "TA 1"=SolvsTA_control_zscore5["TA 1",],
                                 "TA 2"=SolvsTA_control_zscore5["TA 2",])
# - Creating the heatmap -
breakList = seq(min(SolvsTA_control_zscore),max(SolvsTA_control_zscore), by=0.2)
#colorz = colorRampPalette(colorz)(length(breakList))
colorz = scico(n=length(breakList), direction = -1, palette = 'roma')
pheatmap(SolvsTA_control_zscore6, breaks=breakList, color = colorz, cluster_rows = FALSE, show_colnames = FALSE, 
         show_rownames = TRUE, treeheight_row = 0, treeheight_col = 0)

# ---- end ----