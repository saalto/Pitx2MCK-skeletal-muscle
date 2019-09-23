# ---- Required packages ----
require(readxl)
# ---- end ----

# ---- Dataset references ----
# --- genelist consist of gene symbol and log fold change ---
# --- folder of excel spreadsheets of the GO annotation from MGI ---
# --- Need to make sure that the CSV file order is in the same order as the MGI spreadsheets ---

# --- Reference to Reactome dataset ---
genelist = read.csv("C:/Users/sarah/Documents/RNASeq_analysis/2019_08_19/GeneLists/reactome.csv")
setwd("C:/Users/sarah/Documents/RNASeq_analysis/2019_07_30/MGI")
outputFile = "C:/Users/sarah/Documents/RNASeq_analysis/AnnotatedGeneList_reactome_13.csv"

# --- Reference to coFactor dataset ---
genelist = read.csv("C:/Users/sarah/Documents/RNASeq_analysis/2019_08_19/GeneLists/coFactor.csv")
setwd("C:/Users/sarah/Documents/RNASeq_analysis/2019_08_06/MGI_02/coFactors_MGI/")
outputFile = "C:/Users/sarah/Documents/RNASeq_analysis/AnnotatedGeneList_coFactors_16.csv"

# --- Reference to Transcription Factor dataset ---
genelist = read.csv("C:/Users/sarah/Documents/RNASeq_analysis/2019_08_19/GeneLists/TranscriptionFactors.csv")
setwd("C:/Users/sarah/Documents/RNASeq_analysis/2019_08_06/MGI_02/transcriptionfactors_MGI/")
outputFile = "C:/Users/sarah/Documents/RNASeq_analysis/AnnotatedGeneList_TFs_16.csv"

# ---- end ----

# ---- Leave this code (unless DeBugging) ----
files = list.files(".")

# --- Clear and reset the lists ---
pasted_molecular_function = c()
pasted_biological_function = c()
pasted_cellular_component = c()
pasted_references = c()

# --- Need to convert the gene list to a data frame and order the genes alphebetically---
genelist = data.frame(genelist)
#genelist <- genelist[order(row.names(genelist)),]

# --- For molecular function, biological process, and cellular component GO aspect ---
# --- Extract and create a Classifications column for MF, BP, and CC ---
for (i in 1: length(files)){
  data = read_excel(files[i], sheet = 1)

  # - ignore any terms whose evidence comes from ND (no biological data available), 
  #   NAS (non-traceable author statement), or IEA (Inferred from electronic annotation)
  # & data$Evidence != "IEA"
  data <- data[data$Evidence != "ND" & data$Evidence != "NAS",]
  
  list1 = c()
  list2 = c()
  list3 = c()
  list4 = c()
  list1_1 = c()

  for (func in 1:length(data$Aspect)){
    if (data$Aspect[func] == "Molecular Function"){
      list1 = c(list1, data$`Classification Term`[func])
      list1_1 = c(list1_1, data$`Reference(s)`)
      } else if (data$Aspect[func] == "Biological Process"){    
      list2 = c(list2, data$`Classification Term`[func])
      list1_1 = c(list1_1, data$`Reference(s)`)
      } else if (data$Aspect[func] == "Cellular Component"){  
      list3 = c(list3, data$`Classification Term`[func])
      list1_1 = c(list1_1, data$`Reference(s)`)
      } 
    }
    
  with_NA1 = unlist(list1)
  no_NA1 = with_NA1[!is.na(with_NA1)]
  no_NA1_1 = Reduce(unique, no_NA1)
  pasted_molecular_function[i] = paste(no_NA1_1, collapse = ", ") 
  
  with_NA2 = unlist(list2)
  no_NA2 = with_NA2[!is.na(with_NA2)]
  no_NA2_1 = Reduce(union, no_NA2)
  pasted_biological_function[i] = paste(no_NA2_1, collapse = ", ")  
  
  with_NA3 = unlist(list3)
  no_NA3 = with_NA3[!is.na(with_NA3)]
  no_NA3_1 = Reduce(union, no_NA3)
  pasted_cellular_component[i] = paste(no_NA3, collapse = ", ") 
  
  with_NA4 = unlist(list1_1)
  no_NA4 = with_NA4[!is.na(with_NA4)]
  no_NA5 = Reduce(union, no_NA4)
  pasted_references[i] = paste(no_NA5, collapse = "; ")

}

genelist$MolecularFunction = pasted_molecular_function
genelist$BiologicalFunction = pasted_biological_function
genelist$CellularComponent = pasted_cellular_component
genelist$References = pasted_references

# --- Write to file
write.csv(genelist, file = outputFile)

# ---- end ----