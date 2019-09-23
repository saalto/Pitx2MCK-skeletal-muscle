library(ggplot2)

level_order <- c("Muscle contraction",
                 "Cardiac conduction",
                 "Striated Muscle Contraction",
                 "Phase 0 - rapid depolarisation",
                 "Ion homeostasis",
                 "Ion channel transport",
                 "Ion transport by P-type ATPases",
                 "Reduction of cytosolic Ca++ levels",
                 "Amino acid transport across the plasma membrane",
                 "Transport of inorganic cations/anions and amino acids/oligopeptides",
                 "Metabolism of carbohydrates",
                 "Glycolysis","Glucose metabolism",
                 "Gluconeogenesis",
                 "Glycogen breakdown (glycogenolysis)",
                 "Glycogen metabolism",
                 "FGFR3 ligand binding and activation","FGFR3c ligand binding and activation",
                 "Depolymerisation of the Nuclear Lamina")

group <- unique(test1$data$Gene)

# group basedd on log-fold change and 6 main categories
test3 <- ggplot2::ggplot(data=test1$data, 
                         mapping = ggplot2::aes(x=factor(Gene, level=group),y=factor(categoryID, level=level_order),fill=foldChange)) + 
  ggplot2::geom_tile() + 
  ggplot2::scale_fill_gradientn(colours = c("blue", "white", "red")) +
  ggplot2::xlab(label="Gene") + 
  ggplot2::ylab(label="Reactome Category") +
  ggplot2::theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,hjust = 1), 
        axis.text.y = ggplot2::element_text(angle= 0, vjust = 0, hjust = 1)) + 
  ggplot2::scale_x_discrete()
test3

# no grouping based on log-fold change and 6 main categories
test4 <- ggplot2::ggplot(data=test1$data, 
                         mapping = ggplot2::aes(x=Gene,y=factor(categoryID, level=level_order),fill=foldChange)) + 
  ggplot2::geom_tile() + 
  ggplot2::scale_fill_gradientn(colours = c("blue", "white", "red")) +
  ggplot2::xlab(label="Gene") + 
  ggplot2::ylab(label="Reactome Category") +
  ggplot2::theme_bw()+
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,hjust = 1), 
                 axis.text.y = ggplot2::element_text(angle= 0, vjust = 0, hjust = 1)) + 
  ggplot2::scale_x_discrete()
test4