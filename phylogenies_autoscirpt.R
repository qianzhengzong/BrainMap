library("refGenome")
library("dplyr")
library("preprocessCore")
library("ape")

GetEXPRESSION <- function(rgd,spd,tb) {
  # adding individual's RPKM value to col end of orthoglous tables 
  # x: Stringtie tab file full path. (../../BM3-15_.tab,...)
  species <- substr(spd,2,nchar(spd))
  sample.name <- strsplit(tb,"_")[[1]][1]
  tab.file.path <- paste0(rgd,spd,"/",tb)
  stringtie.table <- read.table(tab.file.path,header = F)
  individual.exp.table <- RemoveDUP(stringtie.table,sample.name)
  join.by <- c(species = "gene")
  names(join.by) <- species
  TmpOrthoglous <- left_join(onetoone.orthoglous,individual.exp.table, by = join.by )
  .GlobalEnv$onetoone.orthoglous <- TmpOrthoglous
}

RemoveDUP <- function(st.tab,sp.name) {
  # Remove duplicated ensembl IDs in gene expression table.
  # st.tab: gene expression table produced by Stringtie; name: Sample name (BM3-15,...)
  st.tab <- st.tab[,c(1,8)]
  colnames(st.tab) <- c("gene",sp.name)
  dup <- st.tab[duplicated(st.tab$gene),]$gene
  if (max(duplicated(st.tab$gene)) != 0) {
    st.tab <- st.tab[-which(st.tab$gene %in% dup),] # keep unique genes
  }
  return(st.tab)
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
onetoone.orthoglous <- orthoglous_IDs

region.dir <- "./tab_thalamus/"
species.dir <- list.files("./tab_thalamus/")
for (spdir in species.dir) {
  tab.files <- list.files(paste0(region.dir,spdir))
  for (tabs in tab.files) {
    GetEXPRESSION(region.dir,spdir,tabs)
  }
}


tha.orthoglous_IDs_noNA <- na.omit(onetoone.orthoglous)
rownames(tha.orthoglous_IDs_noNA) <- tha.orthoglous_IDs_noNA[,1]
tha.ExpMatrix <- tha.orthoglous_IDs_noNA[,-c(1:4)]
tha.SampleID <- c("HB44", "HA44", "HC44", "HD44", "MA17", "MB17", "MC17A", "CHB17", "CHC17", "CHD17", "BA17", "BB17" ,"BC17")
colnames(tha.ExpMatrix) <- tha.SampleID

tha.exp.norm <- normalize.quantiles(as.matrix(tha.ExpMatrix))
colnames(tha.exp.norm)<- tha.SampleID
tha.ExpLogNorm <- log2(tha.exp.norm+1)
tha.SpearmanLogNorm <-1- cor(tha.ExpLogNorm,method = "spearman")

tha.Nj_SpLogNorm <- nj(tha.SpearmanLogNorm)

plot.phylo(tha.Nj_SpLogNorm,type = "unrooted",main = "Thalamus NJ tree",use.edge.length = T,node.pos = 1,show.tip.label = T,
           edge.color = "black", cex = 1,no.margin = T)
nodelabels(cex=0.8)
tiplabels(adj = c(1,0), cex=0.8)
edgelabels(cex =0.8)
tha.edglen <- tha.Nj_SpLogNorm$edge.length
tha.edglen <- as.character(round(tha.edglen,4))
edgelabels(text = tha.edglen, col="black",adj = c(1.6,0),cex = 0.8, frame = "none", bg = NULL, pch = 0.001)

library("ouch")
tha.outree <- ape2ouch(tha.Nj_SpLogNorm,scale = TRUE, branch.lengths = tha.Nj_SpLogNorm$edge.length)
