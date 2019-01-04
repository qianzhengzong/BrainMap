library("refGenome")
library("dplyr")


## define fuctions
FunGetExpression <- function(x) {
  # adding individual's RPKM value to col end of orthoglous tables 
  # x: Stringtie tab file under "./tab_ba9/" path. (BM3-15_.tab,...)
  SampleName <- strsplit(x,"_")[[1]][1]
  FilePath <- paste0("./tab_ba9/",x)
  StringtieTable <- read.table(FilePath,header = F)
  IndividualExpTable <- FunRemoveDup(StringtieTable,SampleName)
  TmpOrthoglous <- left_join(orthoglous_IDs,IndividualExpTable, by = c("bonobo"= "gene"))
  .GlobalEnv$orthoglous_IDs <- TmpOrthoglous
}
FunRemoveDup <- function(tab,name) {
  # Remove duplicated ensembl IDs in gene expression table.
  # tab: gene expression table produced by Stringtie; name: Sample name (BM3-15,...)
  tab <- tab[,c(1,8)]
  colnames(tab) <- c("gene",name)
  dup <- tab[duplicated(tab$gene),]$gene
  if (max(duplicated(tab$gene)) != 0) {
    tab <- tab[-which(tab$gene %in% dup),] # keep unique genes
  }
  return(tab)
}

## we get four species 1:1 orthologous genes ensembl list.
orthogous <- read.csv("/home/qzz/Documents/brain_map_analysis/5r_analysis/4_species_orthogolous_id_type_Unique.txt",header = T)
gene_IDs_list <- grep("ene.stable|type",as.vector(colnames(orthogous)),value = T)
gene_IDs <- orthogous[,gene_IDs_list]
gene_IDs <- gene_IDs[gene_IDs[,3] == "ortholog_one2one",]
gene_IDs <- gene_IDs[gene_IDs[,5] == "ortholog_one2one",]
gene_IDs <- gene_IDs[gene_IDs[,7] == "ortholog_one2one",]
orthoglous_IDs <- gene_IDs[,c(1,2,4,6)]
colnames(orthoglous_IDs) <- c("human","macaque","chimpanzee","bonobo")


## combine individuals' expression value(RPKM) to end column of "orthoglous_IDs"to 
## establish expression matrix
AllFiles<- list.files("./tab_ba9/")
lapply(AllFiles, FUN = FunGetExpression)

## remove orthoglous IDs, "ExpMatrix" contains only expression vales
orthoglous_IDs_noNA <- na.omit(orthoglous_IDs)
rownames(orthoglous_IDs_noNA) <- orthoglous_IDs_noNA[,1]
ExpMatrix <- orthoglous_IDs_noNA[,-c(1:4)]
SampleID <- c("HD3", "HC3", "HA3", "HB3", "MB2", "MA2", "MC2", "CHB2", "CHC2", "CHD2", "BA2", "BB2" ,"BC2")
colnames(ExpMatrix) <- SampleID

##below only keep useful variables
keeplist <- c("orthoglous_IDs","SampleTag","keeplist","orthoglous_IDs_noNA","ExpMatrix")
toberemove <- ls()[-which(ls() %in% keeplist)]
# rm(list = toberemove)

## normalization 
library("preprocessCore")
# tmpLog <- log2(tmp+1)
# tmpNormLog <- normalize.quantiles(tmpLog)

par(mfrow=c(4,4))
for (i in c(1:13)) {
  hist(tmpLogNort[,i])
}

exp.norm <- normalize.quantiles(as.matrix(ExpMatrix))


colnames(exp.norm)<- SampleID
ExpLogNorm <- log2(exp.norm+1)
SpearmanLogNorm <-1- cor(ExpLogNorm,method = "spearman")

library("ape")
Nj_SpLogNorm <- nj(SpearmanLogNorm)

edglen <- Nj_SpLogNorm$edge.length
edglen <- as.character(round(edglen,2))
edgelabels(text = edglen, col="black", pch = 0.01,font=1)


library("ouch")
outree <- ape2ouch(Nj_SpLogNorm,scale = TRUE, branch.lengths = Nj_SpLogNorm$edge.length)




