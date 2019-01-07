library("dplyr")

#  TmpOrthoglous <- left_join(orthoglous_IDs,IndividualExpTable, by = c("human"= "gene"))
## predefined varialbles ----
species <- c("human","macaque","chimpanzee","bonobo")
species.list <- as.list(species)
names(species.list) <- species

## we get four species 1:1 orthologous genes ensembl list "orthoglous_IDs".----
orthogous <- read.csv("/home/qzz/Documents/brain_map_analysis/data/4_species_orthogolous_id_type_Unique.txt",header = T)
gene_IDs_list <- grep("ene.stable|type",as.vector(colnames(orthogous)),value = T)
gene_IDs <- orthogous[,gene_IDs_list]
gene_IDs <- gene_IDs[gene_IDs[,3] == "ortholog_one2one",]
gene_IDs <- gene_IDs[gene_IDs[,5] == "ortholog_one2one",]
gene_IDs <- gene_IDs[gene_IDs[,7] == "ortholog_one2one",]
orthoglous_IDs <- gene_IDs[,c(1,2,4,6)]
colnames(orthoglous_IDs) <- c("human","macaque","chimpanzee","bonobo")
write.csv(orthoglous_IDs,file="/home/qzz/Documents/brain_map_analysis/data/orthoglous_IDs.csv",row.names =F)
orthoglous_IDs <- read.csv("/home/qzz/Documents/brain_map_analysis/data/orthoglous_IDs.csv",header = T)

## all libraries and brain regions table "allLibs". ----
allLibs <- read.csv("/home/qzz/Documents/brain_map_analysis/sample_tables/all_lib.csv",header = T)
rownames(allLibs) <- paste0(allLibs$Seq.name,"-",allLibs$Index)
all.species.tabs <- lapply(species.list, function(Fspecies){
  species.table <- allLibs[allLibs$species == Fspecies,]
  return(species.table)
})

##get a list contains all stringtie data.frame. list names are "BM1-10","BM1-11"... ----
tabDir <- "/home/qzz/Documents/brain_map_analysis/data/tab/"
#tabDir <- "/home/qzz/Documents/brain_map_analysis/5r_analysis/tab_ba9/tab2_ba9/human/"
allTabs <- list.files(tabDir)
tabsFullPath <- as.list(paste0(tabDir,allTabs))
tabsNames <- sapply(tabsFullPath,FunGetTabName <- function(Ftabfullpath){
  ##Ffullpath:full path of tab file:"/full/path/BM1-15_.tab"
  ##return:"BM1-15"
  splitedfullpath <- strsplit(Ftabfullpath,"/")[[1]]
  tabname <- splitedfullpath[length(splitedfullpath)]
  rmsuffix <- strsplit(tabname,"_")[[1]][1]
  return(rmsuffix)
})
names(tabsFullPath) <- tabsNames

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
  noduptab <- Ftab
  if (max(duplicated(Ftab$Gene_ID)) != 0) {
    noduptab <- Ftab[-which(Ftab$Gene_ID %in% dup),] # keep unique genes
  }
  rownames(noduptab) <- noduptab[,1]
  return(noduptab)
}
stie.tabs <- lapply(tabsFullPath, FunReadSTieTabs)

FunExtractExp <- function(FrefineTab,FnormMethod = "TPM", Fid) {
  extras <- c("Gene_ID",FnormMethod)
  extractedtab <- FrefineTab[extras]
  #extracted.vector <- as.vector(extractedtab)
  colnames(extractedtab) <- c("Gene_ID",Fid)
  return(extractedtab)
}
tpm.dataframe <- Map(FunExtractExp,FrefineTab=stie.tabs, FnormMethod = "TPM",Fid = names(stie.tabs))

tpm.dataframe.species <- lapply(species.list, function(Fspecies) {
  tab <- all.species.tabs[[Fspecies]]
  bm.list <- rownames(tab)
  extracted.tabs <- tpm.dataframe[bm.list]
  return(extracted.tabs)
})## extract each species "TPM" expression vector

## remove duplicated gene ids in "orthoglous_IDs"----
### find duplicated gene ids in Stringtie output tables,"dup.ids"
FunFindDup <- function(Ftabfullpath) {
  tab <- read.table(Ftabfullpath,header = F)
  stie.cols <- c("Gene_ID","Gene_name","Reference","Strand","Start","End","Coverage","FPKM","TPM")
  colnames(tab) <- stie.cols
  dup.lines <- data.frame()
  dup <- tab[duplicated(tab$Gene_ID),]$Gene_ID
  if (max(duplicated(tab$Gene_ID)) != 0) {
    dup.lines <- tab[which(tab$Gene_ID %in% dup),] ## extract duplicated lines
  }
  return(dup.lines)
}
dup.list <- lapply(tabsFullPath,FunFindDup)
dup.list.species <- lapply(species.list,function(Fspecies) {
  tab <- all.species.tabs[[Fspecies]]
  bm.list <- rownames(tab)
  extracted.tabs <- dup.list[bm.list]
  return(extracted.tabs)
})
dup.list.species.fourlist <- lapply(species.list,function(Fspecies) {
  dlist <- dup.list.species[[Fspecies]]
  df <- data.frame()
  for (i in c(1:length(dlist))) {
    df <- rbind(df,dlist[[i]])
  }
  return(df)
})
dup.ids <- lapply(species.list,function(Fspecies) {
  df <- dup.list.species.fourlist[[Fspecies]]
  dup.vector <- df[!duplicated(df$Gene_ID),]$Gene_ID
  dup.vector <- as.vector(dup.vector)
  return(dup.vector)
})
# rm dup.ids in "orthoglous_IDs"
orth.dup.rows <- lapply(species.list, function(Fspecies) {
  dup <- dup.ids[[Fspecies]]
  rows <- which(orthoglous_IDs[,Fspecies] %in% dup)
  return(rows)
})
orth.dup.vector <- as.vector(unlist(orth.dup.rows))
orthoglous.ids.nodup <- orthoglous_IDs[-orth.dup.vector,]
rownames(orthoglous.ids.nodup) <- c(1:nrow(orthoglous.ids.nodup))

## get species specific expression matrix, "exp.matrix" ----
FunGetSpeciesExpMatrix <- function(Fspecies) {
  exp.species.list <- tpm.dataframe.species[[Fspecies]]
  othgus.species <- as.character(orthoglous.ids.nodup[,Fspecies])
  othgus.species.tbl <- tbl_df(othgus.species)
  colnames(othgus.species.tbl) <- "Gene_ID"
  for (i in 1:length(exp.species.list)) {
    ex.tbl <- tbl_df(exp.species.list[[i]])
    ex.tbl[,1] <- as.character(ex.tbl$Gene_ID) 
    othgus.species.tbl <- left_join(othgus.species.tbl,ex.tbl,by = "Gene_ID")
  }
  return(othgus.species.tbl)
}

tpm.species.matrix<- lapply(species.list,FunGetSpeciesExpMatrix)

oth.tmp <- tbl_df(orthoglous.ids.nodup)
for (i in species) {
  sp <- i
  exp.tmp <- tpm.species.matrix[[i]]
  gid <- "Gene_ID"
  oth.tmp <- inner_join(oth.tmp,exp.tmp,by=setNames(gid,sp))
}
exp.matrix.tmp <- as.data.frame(oth.tmp)
exp.matrix <- exp.matrix.tmp[,5:ncol(exp.matrix.tmp)]


# analysis----


















## some test scripts----

## clear Environment----  
keep.list <- c("allLibs", "species","species.list","orthoglous.ids.nodup", "stie.tabs","tpm.dataframe.species", 
               "dup.list","all.species.tabs","tpm.species.matrix","exp.matrix")
rm.list <- ls()[-which(ls() %in% keep.list)]
rm(list = rm.list)

