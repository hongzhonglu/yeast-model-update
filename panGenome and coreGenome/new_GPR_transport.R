#Enzyme database preparation for whole uniprot database

#library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(readr)



## paramters to filter the mapping results
pident0 <- 45
bitscore0 <- 100


###match
getMultipleReactionFormula <- function(description, reaction_ko, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length( which (reaction_ko %in%  ko[i]))){
      index <- which (reaction_ko %in%  ko[i])
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{
      
      result[i] <- NA
    }
  }
  return(result)
}
splitAndCombine <- function(gene, rxn,sep0) { ##one rxn has several genes, this function was used to splite the genes
  
  gene <- str_split(gene, sep0)
  tt<- length(gene)
  gene0 <- list()
  for (i in 1:tt){
    gene0[[i]] <- paste(rxn[i], gene[[i]], sep = "@@@")
    
  }
  
  gene1 <- unique(unlist(gene0))
  gene2 <- str_split(gene1, "@@@" )
  rxnGene <- data.frame(v1=character(length(gene2)),stringsAsFactors = FALSE)
  tt1 <- length(gene2)
  for (j in 1:tt1){
    rxnGene$v1[j] <- gene2[[j]][2]
    rxnGene$v2[j] <- gene2[[j]][1]
  }
  
  return(rxnGene)
}

##parse the transport protein annotation database
##find transport reactions from http://www.tcdb.org/public/tcdb
transport_enzyme <- scan("Blast results/transport_protein_annotation.txt", sep = "\n", what = "complex")
index0 <- which(str_detect(transport_enzyme,">") ==TRUE)
ss <- vector()
for (i in seq(length(index0))){
  ss[i] <-transport_enzyme[index0[i]]
}


#the relation between gene-protein-reactions were stored in ss
ss <- str_replace_all(ss,">","")
ss0 <- str_split(ss, " ")
ss1 <- data.frame(ID=vector(length=16540), stringsAsFactors = FALSE)

for (i in seq(length(ss0))){
  ss1$ID[i] <- ss0[[i]][1]
  ss1$Annotation[i]  <-  paste(ss0[[i]][-1],collapse = " ")
}
head(ss1)

#further get the ec for each transport annotation
ss2 <- ss1 %>% separate(.,ID, into = c('GNL','Source','proteinID','ec'), sep = "\\|")
ss2$ec <- str_trim(ss2$ec, side = "both")




#input the blast results for pangenome and s288 genome
colname0 <- c("panID",  "ID",    "pident", "length",   "mismatch", "gapopen",  "qstart",  "qend",  "sstart",  "send",  "evalues" , "bitscore")
panGenome <- read.table('Blast results/pangenomeBlastTCDB.txt',header = FALSE, stringsAsFactors = FALSE)
colnames(panGenome) <- colname0
panGenome <- filter(panGenome, pident >= pident0 & bitscore >= bitscore0)
panGenome0 <- panGenome %>% separate(.,ID, into = c('GNL','Source','proteinID','ec'), sep = "\\|")
panGenome$ec <- panGenome0$ec

s288Genome <-read.table('Blast results/s288genomeBlastTCDB.txt',header = FALSE, stringsAsFactors = FALSE)
colnames(s288Genome) <- colname0
s288Genome <- filter(s288Genome, pident >= pident0 & bitscore >= bitscore0)
s288Genome0 <- s288Genome %>% separate(.,ID, into = c('GNL','Source','proteinID','ec'), sep = "\\|")
s288Genome$ec <- s288Genome0$ec


#get the transport based on ec defined in TCDB dabase
newID <- setdiff(panGenome$ec, s288Genome$ec)
index1 <- which(panGenome$ec %in% newID == TRUE)
newsInPan <- panGenome[index1,]
newsInPan$annotation <- getMultipleReactionFormula(ss1$Annotation,ss1$ID,newsInPan$ID)
write.table(newsInPan,"new transport rxn from TCDB.txt", sep = "\t", row.names = FALSE)


#change gene-transport into GPR format
outGPR <- function(rxnID,gene,formula){
  #input:
  #rxnID <- newsInPan$ec
  #gene <- newsInPan$panID
  #formula <- newsInPan$annotation
  #output:
  #result, dataframe to summarize the gene-protein-reaction relation 
  geneAnnotation <- data.frame(gene0=gene,rxn0=rxnID,formula0=formula,stringsAsFactors = FALSE)
  result <- data.frame(reactionID = unique(rxn), stringsAsFactors = FALSE)
  #using getMultipleReactionFormula
  result$gene <- getMultipleReactionFormula(geneAnnotation$gene0, geneAnnotation$rxn0, result$reaction)
  result$formula <- getMultipleReactionFormula(geneAnnotation$formula0, geneAnnotation$rxn0, result$reaction)
  return(result)
  }

newTransporot <-outGPR(newsInPan$ec,newsInPan$panID,newsInPan$annotation)










#transport annotation for s288c genome annotation and analysis
#the s288c genome annotation is quite mature, its annotation has been the reference for a lot of other genomes
#in the present annotation, s288c has the annotation in the TCDB database, while using the blast analysis, we can find
#a lot of hit even we choose the a higher cut-off.
#how to handle this issue?
#the followed analysis showed that if the pidenty is set at the maximum value, then the target can be mapped onto the
#gene sequence of s288c itself.
s288Genome$annotation <- getMultipleReactionFormula(ss1$Annotation,ss1$ID,s288Genome$ID)
s288Genome_pi100 <- filter(s288Genome,pident == 100)

uniprt_tcdb2 <- read_excel("uniprt_tcdb2.xlsx")

#establish one to one relation between gene and tcdb id
sce_tcdb <- splitAndCombine(uniprt_tcdb2$tcdbID,uniprt_tcdb2$gene_name,sep0 = ";")
sce_tcdb <- filter(sce_tcdb, str_detect(sce_tcdb$v1, "\\."))
sce_tcdb$v1 <- str_trim(sce_tcdb$v1, side = "both")
sce_tcdb$type <- 'from_s288c'

#intersection between the blast results and the present annotation
common <- intersect(s288Genome_pi100$panID,sce_tcdb$v2)


#based on the above analysis, the tansport analysis pipeline will be redesigned
#we will focuse on pan-genome non-reference ORF analysis
#the best hit can be divided into three level
#level 0: if the gene has ortholog. It can be merged based on gene of s288c
#level 1: if the hit above the cut off belong to s288c. It can be merged based on gene of s288c
#level 2: if the hit belong to the existed reaction. It can merged based on reactions
#level 3: if the hit belong to the new reaction. It can be regarded as the new GPRs.
panID_with_ortholog <- read_excel("panID_with_ortholog.xlsx")
panID_with_ortholog$type <- 'Ortholog'
nonORF <- panGenome0[1:736,]
nonORF$gene <- getMultipleReactionFormula(uniprt_tcdb2$gene_name,uniprt_tcdb2$Entry,nonORF$proteinID)
nonORF$type <- getMultipleReactionFormula(panID_with_ortholog$type,panID_with_ortholog$panID,nonORF$panID)

#level 0
nonORF_ortholog <- filter(nonORF, type=='Ortholog')


#level 1
nonORF_with288 <- filter(nonORF, !is.na(nonORF$gene) | type =='Ortholog')
nonORF_without288 <- filter(nonORF, is.na(nonORF$gene) & is.na(nonORF$type))
nonORF_without288$ec <- str_trim(nonORF_without288$ec, side = "both")


#level 2 compared with the ec id compared the new ORFs with all the other ORFs
nonORF_without288$source0 <- getMultipleReactionFormula(sce_tcdb$type,sce_tcdb$v1,nonORF_without288$ec)
nonORF_without288$annotation <- getMultipleReactionFormula(ss2$Annotation,ss2$ec,nonORF_without288$ec)
print(length(unique(nonORF_without288$panID)))












