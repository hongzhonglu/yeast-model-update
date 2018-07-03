library(tidyverse)
library(stringr)
library(readxl)

getSingleReactionFormula <- function(description, reaction_ko, ko) {###description can be any charater of metabolite
  index <- vector()
  result <- vector()
  tt <- vector()
  for (i in 1:length(ko)){
    if(length(match(ko[i],reaction_ko))){
      index <- match(ko[i],reaction_ko)
      tt <- description[index]
      result[i] <- paste0(tt, collapse = ";")
    } else{
      
      result[i] <- NA
    }
  }
  return(result)
}
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
AutoUpdate <- function(description1, para1, description2,  para2){ 
  # using the description1 in para1 to update the description2 in para2
  s1 <- list()
  p <- vector()
  description <- vector()
  nn <- length(para2)
  for (i in 1:nn){
    s1[[i]] <- which(para1 %in% para2[i] ==TRUE)
    p[i]<- s1[[i]][1]
  }
  
  for(i in 1:nn){
    if(!is.na(p[i])){
      description[i] <- description1[p[i]]
    } else{
      description[i] <- description2[i]
    }
  }
  return(description)
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

#input the data
gene_yeastGEM <- read_excel("gene_list_yeastGEM.xlsx", sheet = "Sheet1")
pangenome <- scan("allORFs_pangenome.fasta", sep = "\n", what = "complex")

#get the gene name from pangenome file
index <- which(str_detect(pangenome,">")==TRUE)
geneORF <- pangenome[index]
geneORF <- str_replace_all(geneORF,">","")
write.table(geneORF,"geneORF in pangenome.txt", row.names = FALSE, sep = "\t")
panGene <- data.frame(gene = geneORF, stringsAsFactors = FALSE)

##obtain the standard gene name which can be found in S288C
panGene_split <- panGene %>% separate(gene, c("gene","end"), sep ="_Num") %>% 
  separate(gene, c("s1","s2","s3","s4"), sep ="-")
panGene$gene_simple <- paste(panGene_split$s2,panGene_split$s3,panGene_split$s4, sep = "-") %>%
  str_replace_all(.,"-NA","")


# core gene analysis
variable_ORF <- read_excel("1011_yeast_data.xlsx", sheet = "Table S3")
variable_ORF$S288c[1712:2855] <-"S288C"
variable_ORF$gene_type <- "Variable"
variable_ORF0 <- select(variable_ORF,`Annotation Name`)
core_ORF0 <- setdiff(panGene$gene, variable_ORF0$`Annotation Name`)
coreGene <- data.frame(gene = core_ORF0, stringsAsFactors = FALSE)
coreGene$gene_type <- "core_gene"

##obtain the standard gene name which can be found in S288C
coreGene_split <- coreGene %>% separate(gene, c("gene","end"), sep ="_Num") %>% 
  separate(gene, c("s1","s2","s3","s4"), sep ="-")
coreGene$gene_simple <- paste(coreGene_split $s2,coreGene_split $s3,coreGene_split $s4, sep = "-") %>%
  str_replace_all(.,"-NA","")



## mapping analysis
## the S288C gene should be mapping onto panGenome
gene_SGD <- read_excel("yeast_gene_annotation_SGD.xlsx", sheet = "Sheet1")
panGene$SGD_ID <- getMultipleReactionFormula(gene_SGD$SGD_ID,gene_SGD$Systematic_name,panGene$gene_simple)
panGene$gene_type <- getMultipleReactionFormula(coreGene$gene_type,coreGene$gene,panGene$gene)

panGene$gene_type <- AutoUpdate(variable_ORF$gene_type,
                                variable_ORF$`Annotation Name`,
                                panGene$gene_type,
                                panGene$gene)


panGene$exist_s288c <- getMultipleReactionFormula(variable_ORF$S288c,variable_ORF$`Annotation Name`, panGene$gene)
for (i in seq(length(panGene$gene))){
  if(!is.na(panGene$SGD_ID[i])){
    panGene$exist_s288c[i] <-"S288C"
  } else{
    panGene$exist_s288c[i] <- "NA"
  }
}

for (i in seq(length(panGene$gene))){
  if(panGene$gene_type[i] =="core_gene" & panGene$exist_s288c[i] =="NA"){
    panGene$exist_s288c[i] <-"S288C"
  } else{
    panGene$exist_s288c[i] <- panGene$exist_s288c[i]
  }
}

##the metabolic gene from yeastGEM mapping onto pangenome
index1 <- which(panGene$gene_simple %in% gene_yeastGEM$geneNames ==TRUE)
panGene$GEM <- "NO"
panGene$GEM[index1] <- "YES"
write.table(panGene,"panGene.txt", row.names = FALSE, sep = "\t")


# gene in SGD database not in pan-genome
gene_SGD_notPan <- setdiff(gene_SGD$Systematic_name,panGene$gene_simple)
index2 <- which(gene_SGD$Systematic_name %in% gene_SGD_notPan ==TRUE)
gene_SGD_notPan0 <- gene_SGD[index2,] 


#collapsed_ORF amalysis
collapsed_ORF_SGD <- read_excel("panGene_for manual check.xlsx", sheet = "collapsed_ORF_SGD")
collapsed_ORF_SGD$gene_simple <- getSingleReactionFormula(panGene$gene_simple,panGene$gene,collapsed_ORF_SGD$`Annotation Name`)
collapsed_ORF_SGD$exist_in_pan <- collapsed_ORF_SGD$gene_simple

ss1 <- str_split(collapsed_ORF_SGD$`S288c collapsed ORFs`,",")
for (i in seq(length(ss1))){
  collapsed_ORF_SGD$not_in_pan[i] <- paste(setdiff(ss1[[i]],collapsed_ORF_SGD$exist_in_pan[i]), collapse = ";")
}


# summarize genes from s288c not in pangenome due to the collapsed ORF
gene_ortholog <- read_excel("panGene_for manual check.xlsx", sheet = "pan gene with two ortholog SGD")
gene_collapsed1 <- gene_ortholog$not_in_pan %>%
  unique()


collapsed_ORF_SGD <- read_excel("panGene_for manual check.xlsx", sheet = "collapsed_ORF_SGD")
gene_collapsed2 <- collapsed_ORF_SGD$not_in_pan %>% str_split(.,";") %>% 
  unlist() %>%
  unique()


gene_collapsed <- c(gene_collapsed2,gene_collapsed1) %>%
  unique()
gene_collapsed <-  gene_collapsed[which(!is.na(gene_collapsed))]



#merge these collapsed gene with these gene from SGD but not in pan genome
index4 <- which(gene_SGD_notPan0$Systematic_name %in% gene_yeastGEM$geneNames ==TRUE)
gene_SGD_notPan0$GEM[index4] <- "YES"




#mapping the collapsed gene onto panGenome
#these gene from two parts: the ortholog gene and paralog gene(collapsed gene)
#first estabolish the relation between collapsed gene and the gene in panGenome

gene_ortholog <- read_excel("panGene_for manual check.xlsx", sheet = "pan gene with two ortholog SGD")
collapsed_ORF_SGD <- read_excel("panGene_for manual check.xlsx", sheet = "collapsed_ORF_SGD")
gene_needMap1 <- select(gene_ortholog, not_in_pan, exist_in_pan) %>%
  filter(.,!is.na(not_in_pan))
gene_needMap1$panID <- getMultipleReactionFormula(panGene$gene,panGene$gene_simple,gene_needMap1$exist_in_pan)

gene_needMap2 <- splitAndCombine(collapsed_ORF_SGD$not_in_pan,collapsed_ORF_SGD$exist_in_pan,sep0 = ";")
colnames(gene_needMap2) <- c("not_in_pan","exist_in_pan")
gene_needMap2$panID <-getMultipleReactionFormula(panGene$gene,panGene$gene_simple,gene_needMap2$exist_in_pan)

gene_needMap <- rbind.data.frame(gene_needMap1, gene_needMap2)
gene_needMap <- gene_needMap[!duplicated(gene_needMap$not_in_pan),]
gene_needMap$collapsed <- "collapse"
gene_needMap$Remark <- "find panID based on not collapsed gene"




## get the mapped panGene id for these gene not found in panGenome based on the results in the article
gene_SGD_notPan0$panID <- getMultipleReactionFormula(gene_needMap$panID,gene_needMap$not_in_pan,gene_SGD_notPan0$Systematic_name)
gene_SGD_notPan0$collapsed_sign <- getMultipleReactionFormula(gene_needMap$collapsed,gene_needMap$not_in_pan,gene_SGD_notPan0$Systematic_name)
gene_SGD_notPan0$Remark <- getMultipleReactionFormula(gene_needMap$Remark,gene_needMap$not_in_pan,gene_SGD_notPan0$Systematic_name)




## get the mapped panGene id for these gene not found in panGenome based on blast analysis results
blastPanBasedSGD <- read_excel("blastPanBasedSGD for the gene in sgd but no panID.xlsx")


#function to obtain the best panID based on bitscore
getPanID_b <- function(geneID, blast0 = blastPanBasedSGD){
  s1 <- geneID
  if(length(which(blast0$geneID %in% s1==TRUE)) !=0 ){
  index0 <- which(blast0$geneID %in% s1==TRUE)
  target <-blast0[index0,]
  index_maxScore <- which.max(target$bitscore)
  s2 <- target$panID[index_maxScore]
  } else{
  s2 <- "Can't find the panID"
  }
  return(s2)
}

getAAlengthPanID_b <- function(geneID, blast0 = blastPanBasedSGD){
  s1 <- geneID
  if(length(which(blast0$geneID %in% s1==TRUE)) !=0 ){
    index0 <- which(blast0$geneID %in% s1==TRUE)
    target <-blast0[index0,]
    index_maxScore <- which.max(target$bitscore)
    s2 <- target$length[index_maxScore]
  } else{
    s2 <- "Can't find the panID"
  }
  return(s2)
}


#function to obtain the best panID based on pidenty
getPanID_p <- function(geneID, blast0 = blastPanBasedSGD){
  s1 <- geneID
  if(length(which(blast0$geneID %in% s1==TRUE)) !=0 ){
    index0 <- which(blast0$geneID %in% s1==TRUE)
    target <-blast0[index0,]
    index_maxScore <- which.max(target$pident)
    s2 <- target$panID[index_maxScore]
  } else{
    s2 <- "Can't find the panID"
  }
  return(s2)
}


getAAlengthPanID_p <- function(geneID, blast0 = blastPanBasedSGD){
  s1 <- geneID
  if(length(which(blast0$geneID %in% s1==TRUE)) !=0 ){
    index0 <- which(blast0$geneID %in% s1==TRUE)
    target <-blast0[index0,]
    index_maxScore <- which.max(target$pident)
    s2 <- target$length[index_maxScore]
  } else{
    s2 <- "Can't find the panID"
  }
  return(s2)
}




for (i in seq(length(gene_SGD_notPan0$SGD_ID))){
  gene_SGD_notPan0$panID_b[i] <- getPanID_b(gene_SGD_notPan0$Systematic_name[i])
  gene_SGD_notPan0$AAlength_b[i] <- getAAlengthPanID_b(gene_SGD_notPan0$Systematic_name[i])
  gene_SGD_notPan0$panID_p[i] <- getPanID_p(gene_SGD_notPan0$Systematic_name[i])
  gene_SGD_notPan0$AAlength_p[i] <- getAAlengthPanID_p(gene_SGD_notPan0$Systematic_name[i])
  }

gene_SGD_notPan0$CHECK <- gene_SGD_notPan0$panID_b == gene_SGD_notPan0$panID_p








for (i in seq(length(gene_SGD_notPan0$SGD_ID))){
  if(is.na(gene_SGD_notPan0$panID[i]) & gene_SGD_notPan0$CHECK =="TRUE"){
    gene_SGD_notPan0$panID[i] <- gene_SGD_notPan0$panID_b[i]
    gene_SGD_notPan0$Remark[i] <- "find panID based blast"  
  } else{
    gene_SGD_notPan0$panID[i] <- gene_SGD_notPan0$panID[i]
    gene_SGD_notPan0$Remark[i] <- gene_SGD_notPan0$Remark[i]
  }
}


write.table(gene_SGD_notPan0,"gene_SGD_notPan0.txt", sep = "\t",row.names=FALSE)




# Mapping all S288C gene onto panGenomes
gene_SGD$panID <- getMultipleReactionFormula(panGene$gene,panGene$gene_simple,gene_SGD$Systematic_name)
gene_SGD$panID <- AutoUpdate(gene_SGD_notPan0$panID,gene_SGD_notPan0$SGD_ID,gene_SGD$panID,gene_SGD$SGD_ID)
gene_SGD$Remark <- getMultipleReactionFormula(gene_SGD_notPan0$Remark,gene_SGD_notPan0$SGD_ID,gene_SGD$SGD_ID)
write.table(gene_SGD,"gene_SGD_with panID.txt", sep = "\t",row.names=FALSE)







# Ortholog  gene analysis
#in panGenome some panGene have ortholog genes in SGD
#however some ortholog can find the panID from panGenome while some can't find the panID
#ortholog:R0040W, R0020W is not exist in S288C
#YAR033W no panID but it has the paralog

Ortholog1 <- read_excel("panGene_for manual check.xlsx", 
                       sheet = "pan gene with one ortholog SGD")
Ortholog1 <- select(Ortholog1,panID, `Ortholog in SGD_2010`)
colnames(Ortholog1) <- c('panID','systematic_name')
Ortholog1$type <- "single ortholog"

Ortholog2 <- read_excel("panGene_for manual check.xlsx", 
              sheet = "pan gene with two ortholog SGD")
Ortholog2 <- select(Ortholog2,panID, `Ortholog in SGD_2010`)
colnames(Ortholog2) <- c('panID','systematic_name')
Ortholog20 <- splitAndCombine(Ortholog2$systematic_name, Ortholog2$panID, sep=",")
colnames(Ortholog20) <- c('systematic_name','panID')
Ortholog20 <- select(Ortholog20, panID,systematic_name)
Ortholog20$type <- "more orthologs"

panID_ortholog <- rbind.data.frame(Ortholog1, Ortholog20)

#get the blast parameters
panID_ortholog$combine <- paste(panID_ortholog$panID,panID_ortholog$systematic_name, sep = ";")
blastPanBasedSGD$combine <- paste(blastPanBasedSGD$panID, blastPanBasedSGD$geneID,sep = ";")
index5 <- which(blastPanBasedSGD$combine %in% panID_ortholog$combine ==TRUE)
panID_ortholog_BLAST <- blastPanBasedSGD[index5,]


# mapping panGenome onto S288C
# first input the manual check data
panGene$gene_allFrom_s288c <- getMultipleReactionFormula(gene_SGD$Systematic_name,gene_SGD$panID,panGene$gene)
  
  


