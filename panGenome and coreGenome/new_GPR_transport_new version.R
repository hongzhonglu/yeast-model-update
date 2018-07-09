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







#part 2
#new way to annote the gene based on TCDB database
#in this way, we firstly mapping non-ORF in the panID onto the s288c genome based on orthology annotation and 
#blast analysis results.



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
#level 1: if the gene is not ortholog. a best hit will be choosed based on higest bitscore
#if the hit belong to the existed reaction. It can merged based on reactions
#if the hit belong to the new reaction. It can be regarded as the new GPRs.

panID_with_ortholog <- read_excel("panID_with_ortholog.xlsx")
nonORF <- panGenome0[1:736,]
nonORF$gene <- getMultipleReactionFormula(uniprt_tcdb2$gene_name,uniprt_tcdb2$Entry,nonORF$proteinID)
nonORF$ortholog <- getMultipleReactionFormula(panID_with_ortholog$`Ortholog in SGD_2010`,panID_with_ortholog$panID,nonORF$panID)

#level 0
#these gene could mapping onto s288c directly based on ortholog relation
nonORF_ortholog <- filter(nonORF,!is.na(nonORF$ortholog))
panID_ortholog <- unique(nonORF_ortholog$panID) 



#level 1
#these gene could mapping onto s288c based on blast analysis
#new question, based on cut-off: pidenty >= 45, bitscore >=100. One panID could be mapped onto several different
#geneid from s288c, so how to choose the best hit??
#maybe in this step, we could choose just one geneID as the hit for each nonORF based on the bitscore values. It can be found if the panID could be
#mapped onto several genes of s288c itself, then these genes could have the similar function.

panID_notOrtholog <- setdiff(nonORF$panID, panID_ortholog)
index2 <- which(nonORF$panID %in% panID_notOrtholog ==TRUE)
nonORF_notOrtholog <- nonORF[index2,]

#for this panID of not ortholog, we can choose the best hit based on bitscore
#design a code to choose the best hit based on bitscore for each panID

getBestHit <- function (gene, blastList){
#gene1 <- '148-augustus_masked-1999-AID_2'
#blastList <- nonORF_notOrtholog
df1 <- blastList[which(blastList$panID %in% gene),]
ss <- which.max(df1$bitscore)
target <- df1$proteinID[ss]
return(target)
}


for (i in seq(length(nonORF_notOrtholog$panID))){
  nonORF_notOrtholog$best_hit[i] <- getBestHit(nonORF_notOrtholog$panID[i],blastList1)

}


nonORF_notOrtholog$best_hit_s288c <- getMultipleReactionFormula(uniprt_tcdb2$gene_name,uniprt_tcdb2$Entry,nonORF_notOrtholog$best_hit)

#choose the gene with the best hit from s288c
nonORF_notOrtholog_with_hit_s288c <- filter(nonORF_notOrtholog,!is.na(best_hit_s288c))
nonORF_notOrtholog_with_hit_s288c <- nonORF_notOrtholog_with_hit_s288c[!duplicated(nonORF_notOrtholog_with_hit_s288c$panID),]
nonORF_s288c <- select(nonORF_notOrtholog_with_hit_s288c,panID, best_hit_s288c)


#choose the gene with the best hit not from s288c
nonORF_notOrtholog_with_hit_else <- filter(nonORF_notOrtholog, is.na(best_hit_s288c))
nonORF_notOrtholog_with_hit_else <- nonORF_notOrtholog_with_hit_else[!duplicated(nonORF_notOrtholog_with_hit_else$panID),]
nonORF_else <- select(nonORF_notOrtholog_with_hit_else,panID, best_hit)




#level 2 compared with the ec id compared the new ORFs with all the other ORFs
nonORF_else$annotation <- getMultipleReactionFormula(ss2$Annotation,ss2$proteinID,nonORF_else$best_hit)
nonORF_else$ec <- getMultipleReactionFormula(ss2$ec,ss2$proteinID,nonORF_else$best_hit)

#find the gene from s288c based on ec number
nonORF_else$gene_288c_fromEC <- getMultipleReactionFormula(sce_tcdb$v2,sce_tcdb$v1,nonORF_else$ec)

nonORF_newTransport <- select(nonORF_else, annotation, ec)
nonORF_newTransport0 <- nonORF_newTransport[!duplicated(nonORF_newTransport$ec),]
nonORF_newTransport0$panID <- getMultipleReactionFormula(nonORF_else$panID,nonORF_else$ec,nonORF_newTransport0$ec)




#results analysis
#save the results
#total nonORF
gene_nonORF_total <- unique(nonORF$panID)
gene_nonORF_ortholog <- unique(nonORF_ortholog$panID)
gene_nonORF_s288c <- unique(nonORF_s288c$panID)
gene_nonORF_else <- unique(nonORF_else$panID)

x <- c('gene_nonORF_total', 'gene_nonORF_ortholog', 'gene_nonORF_s288c', 'gene_nonORF_else')
y <- c(length(gene_nonORF_total), length(gene_nonORF_ortholog), length(gene_nonORF_s288c), length(gene_nonORF_else))
resulst <- data.frame(type=x, num=y, stringsAsFactors = FALSE)


library(ggplot2)
ggplot(data=resulst, aes(x = reorder(type, -num), y = num)) + 
  geom_bar(stat = "identity") +
  xlab("Gene_type") + ylab("Number of nonORF gene based on TCDB database")



write.table(nonORF_ortholog, "Result_final/nonORF_ortholog_TCDB.txt", row.names = FALSE, sep = "\t") # the orf with ortholog
write.table(nonORF_s288c, "Result_final/nonORF_s288c_TCDB.txt", row.names = FALSE, sep = "\t") # the orf with best hit from s288c
write.table(nonORF_newTransport0, "Result_final/nonORF_newTransport0.txt", row.names = FALSE, sep = "\t") # the orf with the hit from other strains



