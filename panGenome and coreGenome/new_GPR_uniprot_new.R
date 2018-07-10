#Enzyme database preparation for whole uniprot database

#library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(readr)


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

#summary the protein and EC number relation from Expasy
protein_ec <- scan("Blast results/uniprot_sport_annotation.txt", sep = "\n", what = "complex")
index0 <- which(str_detect(protein_ec,"//") ==TRUE)
index0 <- index0[-1]
ss <- list()
for (i in seq(length(index0)-1)){
  ss[[i]] <- protein_ec[(index0[i]+1):(index0[i+1]-1)]
}
#note
#the relation between gene-protein-reactions were stored in ss
#with the uniprot id of protein, we can get the related ec and reaction


#function to get the ec number based on uniprot protein iD
#Q07121 "ID   1.4.3.21"
getEC <- function(id,s0=ss){
  require(stringr)
  exist <- vector()
  for (i in seq(length(s0))){
    exist[[i]] <- any(str_detect(s0[[i]],id) == TRUE)
    
  }
  ec <- vector()
  if(length(which(exist==TRUE))){
    target <- which(exist==TRUE)
    #target <- c(961,1000)
    ec_annotation <- s0[target]
    
    for (i in seq(length(ec_annotation))){
      ec[i] <- ec_annotation[[i]][1]
    }
    
    ec_combine <- paste0(ec,collapse = " ", sep = ";")
  } else{
    ec_combine <- NA
  }
  return(ec_combine)
}

getEC(id="Q07121") # just an example

# these two function: getRxn and getRXNstring could be used to get the rxn information from uniprot database
getRxn <- function(id,s0=ss){
  require(stringr)
  #id="Q9I492"
  exist <- vector()
  for (i in seq(length(s0))){
    exist[[i]] <- any(str_detect(s0[[i]],id) == TRUE)
    
  }
  rxn <- vector()
  if(length(which(exist==TRUE))){
    target <- which(exist==TRUE)
    #target <- c(961,1000)
    ec_annotation <- s0[target]
   
    for (i in seq(length(ec_annotation))){
      rxn[i] <- getRXNstring(ec_annotation[[i]])
    }
    
    rxn_combine <- paste0(rxn,collapse = " ", sep = ";")
  } else{
    rxn_combine <- NA
  }
  
  return(rxn_combine)
}

getRXNstring <- function(ss){
  tt <- unlist(ss)
  index0 <- which(str_detect(tt,"^CA") ==TRUE)
  if(length(index0)){
    mm <- tt[index0]
    mm <- str_replace_all(mm, "^CA", "")
    rxn_string <- paste(mm, collapse = ' ')
  } else{
    rxn_string <- ''
  }
  
  return(rxn_string)
}

getRxn(id="Q07121") # just an example

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






#input and filter the data
## paramters to filter the mapping results
pident0 <- 45
bitscore0 <- 100


#Blast_with_uniprot <- read_excel("Blast results/Blast with uniprot.xlsx")

ss0 <- c("geneID","uniprotID", "pident", "length", "mismatch","gapopen","qstart", "qend", "sstart", "send", "evalues", "bitscore" )

panGenome <- read.table('Blast results/pangenomeBlastUniprot.txt',header = FALSE, stringsAsFactors = FALSE)
colnames(panGenome) <- ss0
panGenome <- filter(panGenome, pident >= pident0 & bitscore >= bitscore0)

s288Genome <-read.table('Blast results/s288genomeBlastUniprot.txt',header = FALSE, stringsAsFactors = FALSE)
colnames(s288Genome) <- ss0
s288Genome <- filter(s288Genome, pident >= pident0 & bitscore >= bitscore0)

newID <- setdiff(panGenome$uniprotID, s288Genome$uniprotID)
index1 <- which(panGenome$uniprotID %in% newID == TRUE)
newsInPan <- panGenome[index1,]

# split the uniprotID for the mapping
panGene0 <- newsInPan %>% separate(uniprotID, into = c('sp','uniprotID','gene_name'), sep = "\\|")

# find the EC number for the choosed panGene0
for (i in 1:length(panGene0$geneID)){
  panGene0$EC[i] <- getEC(panGene0$uniprotID[i])
  panGene0$rxn[i] <- getRxn(panGene0$uniprotID[i])
  }



#find reaction based on ec number from metnet database
#establish the gene-ec relation
panGene0 <- panGene0[!is.na(panGene0$EC),]

panGene0_ec <- splitAndCombine(panGene0$EC,panGene0$geneID,sep0 = ";")

colnames(panGene0_ec) <- c('EC','panID')
panGene0_ec <- select(panGene0_ec, panID, EC)
panGene0_ec <- filter(panGene0_ec, str_detect(panGene0_ec$EC,"ID"))
panGene0_ec$EC <- str_replace_all(panGene0_ec$EC, "ID ","")


#input the original ec in s288c
uniprot_s288c <- read_excel("uniprot_s288c.xlsx")
ec_s288c <- select(uniprot_s288c, Entry,`EC number`)
colnames(ec_s288c) <- c('uniprotID','EC')
ec_s288c0 <- splitAndCombine(ec_s288c$EC,ec_s288c$uniprotID,sep0 = ";")
colnames(ec_s288c0) <- c('EC','uniprotID')
ec_s288c0 <- select(ec_s288c0,uniprotID,EC) 
ec_s288c0 <- filter(ec_s288c0,str_detect(ec_s288c0$EC,"\\."))

panGene0_ec$EC <- str_trim(panGene0_ec$EC, side = "both")
ec_s288c0$EC <- str_trim(ec_s288c0$EC, side = "both")
newEC1 <- setdiff(panGene0_ec$EC,ec_s288c0$EC)
panGene0_newEC <- panGene0_ec[which(panGene0_ec$EC %in% newEC1),]

unique(panGene0_newEC$panID)
unique(panGene0_newEC$EC)


getDetailRxnFromMetnet <-function (dataframe0){
  #input
  #dataframe0: column 1, panID; column 2, EC; 
  
  #get the gene-protein-reaction
  metanetx_rxn <- read_tsv("reac_prop_metanetx.tsv")
  metanetx_rxn0 <- select(metanetx_rxn,MNX_ID,EC)
  metanetx_rxn1 <- splitAndCombine(metanetx_rxn0$EC,metanetx_rxn0$MNX_ID,sep0 = ";")
  metanetx_rxn1$v1 <- str_trim(metanetx_rxn1$v1, side = "both")
  
  #dataframe0 <- panGene0_newEC #test example
  dataframe0$MNXID <- getMultipleReactionFormula(metanetx_rxn1$v2,metanetx_rxn1$v1, dataframe0$EC)
  dataframe0 <-  dataframe0[which(!is.na( dataframe0$MNXID)),]
  panGene_rxn <- splitAndCombine( dataframe0$MNXID, dataframe0$panID,sep0 = ";")
  panGene_rxn$EC <- getMultipleReactionFormula( dataframe0$EC, dataframe0$panID,panGene_rxn$v2)
  panGene_rxn0 <- data.frame(MNXID=unique(panGene_rxn$v1))
  panGene_rxn0$panID <- getMultipleReactionFormula(panGene_rxn$v2,panGene_rxn$v1,panGene_rxn0$MNXID)
  panGene_rxn0$EC <- getMultipleReactionFormula(panGene_rxn$EC,panGene_rxn$v1,panGene_rxn0$MNXID)
  ss_ec <- str_split(panGene_rxn0$EC, ";")
  for (i in seq(length(ss_ec))){
    ss_ec[[i]] <- unique(ss_ec[[i]])
    panGene_rxn0$EC[i] <- paste(ss_ec[[i]],collapse = ";")
  }
  panGene_rxn0$Equation <- getMultipleReactionFormula(metanetx_rxn$Equation,metanetx_rxn$MNX_ID,panGene_rxn0$MNXID)
  panGene_rxn0$Description <- getMultipleReactionFormula(metanetx_rxn$Description,metanetx_rxn$MNX_ID,panGene_rxn0$MNXID)
  
  return(panGene_rxn0)
}


panGene_rxn0 <- getDetailRxnFromMetnet(panGene0_newEC)




#save1
outdir <- paste('Blast results/new rxn from pan based on uniprot',pident0,bitscore0,sep = "_")
outdir0 <- paste(outdir,'.txt', sep = "")
write.table(panGene_rxn0, outdir0, row.names = FALSE, sep = "\t")


#save2
#write.table(panGene_rxn0,'Blast results/new rxn from pan based on uniprot.txt', row.names = FALSE, sep = "\t")
#write.table(panGene0_newEC,'Blast results/new ec from pan based on uniprot.txt', row.names = FALSE, sep = "\t")










##method 2
#based on the above analysis, the tansport analysis pipeline will be redesigned
#we will focuse on pan-genome non-reference ORF analysis
#the best hit can be divided into three level
#level 0: if the gene has ortholog. It can be merged based on gene of s288c
#level 1: if the hit above the cut off belong to s288c. It can be merged based on gene of s288c
#level 2: if the hit belong to the existed reaction. It can merged based on reactions
#level 3: if the hit belong to the new reaction. It can be regarded as the new GPRs

##s288 genome annotation could get from SGD directly
##then the compartion could be easy
##non-ORF will analyse independently
panGene_for_manual_check <- read_excel("panGene_for manual check.xlsx")
nonORF_ID <- panGene_for_manual_check$gene[1:1715]
index_nonORF <- which(panGenome$geneID %in% nonORF_ID ==TRUE)
nonORF <- panGenome[index_nonORF,]
nonORF <- nonORF %>% separate(.,uniprotID, into = c('sp','uniprotID0','source'), sep = "\\|")


panID_with_ortholog <- read_excel("panID_with_ortholog.xlsx")
panID_with_ortholog$type <- 'Ortholog'
uniprt_tcdb <- read_excel("uniprt_tcdb2.xlsx")
uniprt_tcdb$type <- 'sce'

nonORF$type <- getMultipleReactionFormula(panID_with_ortholog$type,panID_with_ortholog$panID,nonORF$geneID)
nonORF$gene <- getMultipleReactionFormula(uniprt_tcdb$gene_name,uniprt_tcdb$Entry,nonORF$uniprotID0)



#level 0
#these gene could mapping onto s288c directly based on ortholog relation
nonORF_ortholog <- filter(nonORF, type=='Ortholog')
panID_ortholog <- unique(nonORF_ortholog$geneID)


#level 1
#these gene could mapping onto s288c based on blast analysis
#new question, based on cut-off: pidenty >= 45, bitscore >=100. One panID could be mapped onto several different
#geneid from s288c, so how to choose the best hit??

panID_nonORF <- unique(nonORF$geneID)
panID_notOrtholog <- setdiff(panID_nonORF, panID_ortholog)
index2 <- which(nonORF$geneID %in% panID_notOrtholog ==TRUE)
nonORF_notOrtholog <- nonORF[index2,]


#for this panID of not ortholog, we can choose the best hit based on bitscore
#design a code to choose the best hit based on bitscore for each panID

getBestHit <- function (gene, blastList){
  #gene <- '4-EC1118_1F14_0045g'
  #blastList <- nonORF_notOrtholog
  df1 <- blastList[which(blastList$geneID %in% gene),]
  ss <- which.max(df1$bitscore)
  target <- df1$uniprotID0[ss]
  return(target)
}


for (i in seq(length(nonORF_notOrtholog$geneID))){
  nonORF_notOrtholog$best_hit[i] <- getBestHit(nonORF_notOrtholog$geneID[i], nonORF_notOrtholog)
}

nonORF_notOrtholog$best_hit_s288c <- getMultipleReactionFormula(uniprt_tcdb$gene_name,uniprt_tcdb$Entry,nonORF_notOrtholog$best_hit)


#choose the gene with the best hit from s288c
nonORF_notOrtholog_with_hit_s288c <- filter(nonORF_notOrtholog,!is.na(best_hit_s288c))
nonORF_notOrtholog_with_hit_s288c <- nonORF_notOrtholog_with_hit_s288c[!duplicated(nonORF_notOrtholog_with_hit_s288c$geneID),]
nonORF_s288c <- select(nonORF_notOrtholog_with_hit_s288c,geneID, best_hit_s288c)
length(unique(nonORF_s288c$geneID))

#choose the gene with the best hit not from s288c
nonORF_notOrtholog_with_hit_else <- filter(nonORF_notOrtholog, is.na(best_hit_s288c))
nonORF_notOrtholog_with_hit_else <- nonORF_notOrtholog_with_hit_else[!duplicated(nonORF_notOrtholog_with_hit_else$geneID),]
nonORF_else <- select(nonORF_notOrtholog_with_hit_else,geneID, best_hit)



#level 2 compared with the ec id 
#compared the new ORFs with all the other ORFs
#these gene can be merged if they have the similar ec number with the s288c genome
#if not, then the new gene-protein-rxn could be found for these panIDs
# find the EC number for the choosed panGene0
for (i in 1:length(nonORF_else$geneID)){
  nonORF_else$EC[i] <- getEC(nonORF_else$best_hit[i])
  nonORF_else$rxn[i] <- getRxn(nonORF_else$best_hit[i])
}

# find rhea reaction
rhea2uniprot <- read_tsv('rhea2uniprot.tsv')
nonORF_else$rxn_rhea <- getMultipleReactionFormula(rhea2uniprot$MASTER_ID,rhea2uniprot$ID,nonORF_else$best_hit)




# find new EC based on previous method
nonORF_else <- nonORF_else[!is.na(nonORF_else$EC),]
nonORF_else_ec <- splitAndCombine(nonORF_else$EC,nonORF_else$geneID,sep0 = ";")

colnames(nonORF_else_ec) <- c('EC','panID')
nonORF_else_ec <- select(nonORF_else_ec, panID, EC)
nonORF_else_ec <- filter(nonORF_else_ec, str_detect(nonORF_else_ec$EC,"ID"))
nonORF_else_ec$EC <- str_replace_all(nonORF_else_ec$EC, "ID ","") %>%
  str_trim(., side = "both")

newEC2 <- setdiff(nonORF_else_ec$EC, ec_s288c0$EC)
index_newEC <- which(nonORF_else_ec$EC %in% newEC2)
panID_withNewEc <- nonORF_else_ec[index_newEC,]

newRxn_for_merge <- getDetailRxnFromMetnet(panID_withNewEc)
panID_withNewEc$rxn <- getMultipleReactionFormula(nonORF_else$rxn,nonORF_else$geneID,panID_withNewEc$panID)
panID_withNewEc$rxn_rhea <- getMultipleReactionFormula(nonORF_else$rxn_rhea,nonORF_else$geneID,panID_withNewEc$panID)


#if not new ec, then we could establish the mapping between this gene with s288c genome
commonEC <- intersect(nonORF_else_ec$EC, ec_s288c0$EC)
index_commonEC<- which(nonORF_else_ec$EC %in% commonEC)
panID_commonEC <- nonORF_else_ec[index_commonEC,]
panID_commonEC$rxn <- getMultipleReactionFormula(nonORF_else$rxn,nonORF_else$geneID,panID_commonEC$panID)
panID_commonEC$rxn_rhea <- getMultipleReactionFormula(nonORF_else$rxn_rhea,nonORF_else$geneID,panID_commonEC$panID)

panID_commonEC$geneS288c <- getMultipleReactionFormula(ec_s288c0$uniprotID,ec_s288c0$EC,panID_commonEC$EC)




#we can further establish mapping between the panID with the gene from s288c if they have the common ec, but in this case we could obtain the
#the gene mapping relation from kegg web service, which will be the key!


#save the new results
write.table(panID_withNewEc,'Result_final/newRxn_for_merge from uniprot database.txt', sep = "\t", row.names = FALSE)
write.table(newRxn_for_merge,'Result_final/newRxn_for_merge from uniprot database_detail.txt', sep = "\t", row.names = FALSE)







