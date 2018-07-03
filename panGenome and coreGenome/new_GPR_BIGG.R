#library(tidyverse)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)


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
splitAndCombine0 <- function(gene, rxn, sep0) { ##one rxn has several genes, this function was used to splite the genes
  
  gene <- str_split(gene, sep0)
  tt<- length(gene)
  gene0 <- list()
  for (i in 1:tt){
    gene0[[i]] <- paste(rxn[i], gene[[i]], sep = "@@@")
    
  }
  
  gene1 <- unlist(gene0)
  gene2 <- str_split(gene1, "@@@" )
  rxnGene <- data.frame(v1=character(length(gene2)),stringsAsFactors = FALSE)
  tt1 <- length(gene2)
  for (j in 1:tt1){
    rxnGene$v1[j] <- gene2[[j]][2]
    rxnGene$v2[j] <- gene2[[j]][1]
  }
  
  return(rxnGene)
}

getwd()
gprs_bigg <- scan("BiGG/bigg_proteins.faa", sep = "\n", what = "complex")

head(gprs_bigg)
index0 <- which(str_detect(gprs_bigg,">") ==TRUE)

ss <- vector()
for (i in seq(length(index0))){
  ss[i] <- gprs_bigg[index0[i]]
}

ss <- str_replace_all(ss,">","")
head(ss)

bigg_gprs <- read_excel("BiGG/bigg_gprs.xlsx")

for (i in seq(length(bigg_gprs$gene))){
  bigg_gprs$protein_id[i] <- str_replace_all(bigg_gprs$gene[i],'G_',paste(bigg_gprs$model[i],".", sep = ""))

}


##merge the latest information iML1515 and Recon3D
gpr_human <- read_excel("BiGG/gpr_human.xlsx", 
                        sheet = "Sheet1")
gpr_ecoli <- read_excel("BiGG/gpr_ecoli.xlsx", 
                        sheet = "iML1515_GP")


    # gpr summary in human
human_protein <- scan("BiGG/human_protein.fasta", sep = "\n", what = "complex")
gpr_human <- read_excel("BiGG/gpr_human.xlsx")

summaryLatestGPR <- function(gprs_fasta, gpr){
   index1 <- which(str_detect(gprs_fasta,">") ==TRUE)
   hm <- vector()
   for (i in seq(length(index1))){
      hm[i] <- gprs_fasta[index1[i]]
   }

   hm <- str_replace_all(hm,">","")
   hm <- str_split(hm, " ")
   hm1 <- vector()
   for (i in seq(length(hm))){
      hm1[i] <- hm[[i]][1]
   }

   hm2 <- data.frame(ID1=hm1, ID=hm1, stringsAsFactors = FALSE)
   hm20 <- hm2 %>% separate(ID, into = c("sp","UniprotID","source"), sep = "\\|")

   hm20$rxn <- getMultipleReactionFormula(gpr$m_reaction,gpr$seq_uniprot,hm20$UniprotID)

   hm20 <- hm20[!is.na(hm20$rxn),]

   hm3<- splitAndCombine0(hm20$rxn,hm20$ID1, sep0 = ";")
   colnames(hm3) <- c('reaction','protein_id')
   return(hm3)
}

gpr_human0 <- summaryLatestGPR(human_protein,gpr_human)


    # gpr summary for e.coli
ecoli_protein <- scan("BiGG/ecoli_protein.fasta", sep = "\n", what = "complex")
gpr_ecoli <- read_excel("BiGG/gpr_ecoli.xlsx")
gpr_ecoli0 <- summaryLatestGPR(ecoli_protein,gpr_ecoli)


#merge the gpr-bigg with the latest e.coli and human model
bigg_gprs0 <- select(bigg_gprs,reaction, protein_id)
bigg_gprs1 <- rbind.data.frame(bigg_gprs0,gpr_human0,gpr_ecoli0) # contains the latest genome annotation from e.coli and human



#input the blast results for pangenome and s288 genome
ss <- c("panID",  "ID",    "pident", "length",   "mismatch", "gapopen",  "qstart",  "qend",  "sstart",  "send",  "evalues" , "bitscore")
panGenome <- read.table('BiGG/pangenomeBlastBiGG.txt',header = FALSE, stringsAsFactors = FALSE)
colnames(panGenome) <- ss
panGenome <- filter(panGenome, pident >= pident0 & bitscore >= bitscore0)

s288Genome <-read.table('BiGG/s288genomeBlastBiGG.txt',header = FALSE, stringsAsFactors = FALSE)
colnames(s288Genome) <- ss
s288Genome <- filter(s288Genome, pident >= pident0 & bitscore >= bitscore0)

newID <- setdiff(panGenome$ID, s288Genome$ID)
index1 <- which(panGenome$ID %in% newID == TRUE)
newsInPan <- panGenome[index1,]


#input the bigg reactions
bigg_models_reactions <- read_excel("BiGG/bigg_models_reactions.xlsx")

#function to get the result
getRxnBlast <- function(blast_result,
                        bigg_rxn=bigg_models_reactions,
                        bigg_gprs=bigg_gprs1){
  # blast result is from diamond result
   
  #Blast_with_bigg0 <- newsInPan
  # input the reaction information
  Blast_with_bigg0$rxn <- getMultipleReactionFormula(bigg_gprs$reaction,bigg_gprs$protein_id,Blast_with_bigg0$ID)

  Blast_with_bigg0 <- filter(Blast_with_bigg0, !is.na(rxn))

  Blast_with_bigg1 <- select(Blast_with_bigg0, panID,ID,rxn)
  colnames(Blast_with_bigg1) <- c('panID','geneID','reactionID')

  # get the reaction based on the reactonID
  rxn_bigg <- splitAndCombine0(Blast_with_bigg1$reactionID,Blast_with_bigg1$panID,sep0 = ";")
  rxn_bigg$v1 <- str_trim(rxn_bigg$v1, side = 'both')
  rxn_bigg$v2 <- str_trim(rxn_bigg$v2, side = 'both')
  
  rxn_bigg$v3 <- paste(rxn_bigg$v1,rxn_bigg$v2)
  rxn_bigg <- rxn_bigg[!duplicated(rxn_bigg$v3),]
  colnames(rxn_bigg) <- c('biggID','panID','MIX')
  rxn_bigg$biggID <- str_replace_all(rxn_bigg$biggID, "R_","")
  
  rxn_bigg$biggID <- str_replace_all(rxn_bigg$biggID,'HM','HMR_')
  rxn_bigg$rxn <- getMultipleReactionFormula(bigg_rxn$reaction_string,bigg_rxn$bigg_id,rxn_bigg$biggID)
  
  return(rxn_bigg)

}

newrxnPan <- getRxnBlast(newsInPan) # it should be noted that due to some differences in the id for the reactions, some rxns should be manually checked.

write.table(newrxnPan,'new rxn from pan.txt', sep = "\t", row.names = FALSE)





