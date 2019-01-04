library("dplyr")

#  TmpOrthoglous <- left_join(orthoglous_IDs,IndividualExpTable, by = c("human"= "gene"))
## predefined varialbles ----
all.species <- c("human","macaque","chimpanzee","bonobo")

## all defined functions. ----

FunReadSTieTabs <- function(Ftabfullpath) {
  ##Ffullpath:full path of tab file:"/full/path/BM1-15_.tab"
  ##return: refined stringtie tables
  tab <- read.table(Ftabfullpath,header = F)
  refinedtab <- FunRefineSTieTab(tab)
  id.keep.list <- c("Gene_ID","FPKM","TPM") 
    #alternatives,c("Gene_ID","Gene_name","Reference","Strand","Start","End","Coverage","FPKM","TPM")
  refinedtab <- refinedtab[id.keep.list]
  return(refinedtab)
}
FunRefineSTieTab <- function(Ftab) {
  # Add colnames and rownames,Remove duplicated ensembl IDs in gene expression table.
  # Ftab: gene expression table produced by Stringtie; name: Sample name (BM3-15,...)
  # return: Refined stringtie tables.
  stie.cols <- c("Gene_ID","Gene_name","Reference","Strand","Start","End","Coverage","FPKM","TPM")
  colnames(Ftab) <- stie.cols
  dup <- Ftab[duplicated(Ftab$Gene_ID),]$Gene_ID
  if (max(duplicated(Ftab$Gene_ID)) != 0) {
    noduptab <- Ftab[-which(Ftab$Gene_ID %in% dup),] # keep unique genes
  }
  rownames(noduptab) <- noduptab[,1]
  return(noduptab)
}
FunGetTabName <- function(Ftabfullpath){
  ##Ffullpath:full path of tab file:"/full/path/BM1-15_.tab"
  ##return:"BM1-15"
  splitedfullpath <- strsplit(Ftabfullpath,"/")[[1]]
  tabname <- splitedfullpath[length(splitedfullpath)]
  rmsuffix <- strsplit(tabname,"_")[[1]][1]
  return(rmsuffix)
}
FunExtractExp <- function(FrefineTab,FnormMethod = "TPM") {
  extras <- c(FnormMethod)
  extractedtab <- FrefineTab[extras]
#  extracted.vector <- as.vector(extractedtab)
  return(extractedtab)
}
FunGetSpeciesExpMatrix <- function() {
  
}

## we get four species 1:1 orthologous genes ensembl list "orthoglous_IDs".----
orthogous <- read.csv("/home/qzz/Documents/brain_map_analysis/5r_analysis/4_species_orthogolous_id_type_Unique.txt",header = T)
gene_IDs_list <- grep("ene.stable|type",as.vector(colnames(orthogous)),value = T)
gene_IDs <- orthogous[,gene_IDs_list]
gene_IDs <- gene_IDs[gene_IDs[,3] == "ortholog_one2one",]
gene_IDs <- gene_IDs[gene_IDs[,5] == "ortholog_one2one",]
gene_IDs <- gene_IDs[gene_IDs[,7] == "ortholog_one2one",]
orthoglous_IDs <- gene_IDs[,c(1,2,4,6)]
colnames(orthoglous_IDs) <- c("human","macaque","chimpanzee","bonobo")

## all libraries and brain regions table "allLibs". ----
allLibs <- read.csv("/home/qzz/Documents/brain_map_analysis/sample_tables/all_lib.csv",header = T)
rownames(allLibs) <- paste0(allLibs$Seq.name,"-",allLibs$Index)
allLibs.species.names <- as.list(all.species)
names(allLibs.species.names) <- all.species
all.species.tabs <- lapply(allLibs.species.names, function(Fspecies){
  species.table <- allLibs[allLibs$species == Fspecies,]
  return(species.table)
})

##get a list contains all stringtie data.frame. list names are "BM1-10","BM1-11"... ----
tabDir <- "/home/qzz/Documents/brain_map_analysis/5r_analysis/tab/"
#tabDir <- "/home/qzz/Documents/brain_map_analysis/5r_analysis/tab_ba9/tab2_ba9/human/"
allTabs <- list.files(tabDir)
tabsFullPath <- as.list(paste0(tabDir,allTabs))
tabsNames <- sapply(tabsFullPath,FunGetTabName)
names(tabsFullPath) <- tabsNames

stie.tabs <- lapply(tabsFullPath, FunReadSTieTabs)
tpm.dataframe <- lapply(stie.tabs, FunExtractExp,"TPM") ## optional "RPKM","Coverage"


# allLibs.species.names <- allLibs.species.names[allLibs.species.names != "macaque"]
tpm.dataframe.species <- lapply(allLibs.species.names, function(Fspecies) {
  tab <- all.species.tabs[[Fspecies]]
  bm.list <- rownames(tab)
  extracted.tabs <- tpm.dataframe[bm.list]
  return(extracted.tabs)
})## extract each species "TPM" expression vector

# tpm.dataframe.orthouglous <- lapply(tpm.dataframe.species,function)

rm.list = ls()[ls() != "stie.tabs"]

