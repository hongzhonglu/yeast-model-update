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


# panGenome annotation based on : https://www.genome.jp/tools/kaas/  
# use BBH (bi-directional best hit) to get the annotation
panGenome_kegg <- read.delim2("panGenome_annotation_egg/panProtein_annotation_kegg.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(panGenome_kegg) <- c('query','ko')

panGenome_kegg0 <- filter(panGenome_kegg, ko != "")



#input the KO annotation from kegg database
#KO_rxn mapping
KO_rxn <- read.delim2("panGenome_annotation_egg/reaction-KO from kegg .txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(KO_rxn) <- c('rxnID','ko')
KO_rxn$ko <- str_replace_all(KO_rxn$ko, "ko:", "")
# rxn formual in kegg
rxn_kegg <- read.delim2("panGenome_annotation_egg/reaction_kegg summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#KO_pathway mapping
KO_pathway <- read.delim2("panGenome_annotation_egg/ko_pathway.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(KO_pathway) <- c('pathway','ko')
KO_pathway$ko <- str_replace_all(KO_pathway$ko, "ko:", "")
KO_pathway <- filter(KO_pathway,str_detect(KO_pathway$pathway,"ko")==FALSE)
#pathway annotation
pathway <- read.delim2("panGenome_annotation_egg/pathway_list_kegg.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(pathway) <- c('pathway','pathway_name')


#get the rxn based on KO
panGenome_kegg0$rxns <- getMultipleReactionFormula(KO_rxn$rxnID,KO_rxn$ko,panGenome_kegg0$ko)
panGenome_newRxn <- filter(panGenome_kegg0,!is.na(rxns))

panGenome_newRxn0 <- splitAndCombine(panGenome_newRxn$rxns, panGenome_newRxn$query, sep0 = ";")
colnames(panGenome_newRxn0) <- c('rxn','query')

panGenome_newRxn0$formula <- getMultipleReactionFormula(rxn_kegg$reaction,rxn_kegg$rxnID,panGenome_newRxn0$rxn)
panGenome_newRxn0$description <- getMultipleReactionFormula(rxn_kegg$name,rxn_kegg$rxnID,panGenome_newRxn0$rxn)
panGenome_newRxn0$ko <- getMultipleReactionFormula(panGenome_kegg0$ko,panGenome_kegg0$query,panGenome_newRxn0$query)
panGenome_newRxn0$pathway <- getMultipleReactionFormula(KO_pathway$pathway,KO_pathway$ko,panGenome_newRxn0$ko)

gene_pathway <- splitAndCombine(panGenome_newRxn0$pathway, panGenome_newRxn0$query, sep0 = ";")
colnames(gene_pathway)  <- c('pathID','gene')
gene_pathway$pathway_name <- getMultipleReactionFormula(pathway$pathway_name,pathway$pathway,gene_pathway$pathID)
#change the pathway into sce
gene_pathway$pathID <- str_replace_all(gene_pathway$pathID,"path:map","sce")
gene_pathway$pathway_name0 <- paste(gene_pathway$pathID, gene_pathway$pathway_name, sep = "  ")
#merge the pathways
panGenome_newRxn0$pathway_name0 <- getMultipleReactionFormula(gene_pathway$pathway_name0,gene_pathway$gene,panGenome_newRxn0$query)



#obtain the new Rxns
newOK <- setdiff(panGenome_newRxn0$rxn[1:501], panGenome_newRxn0$rxn[502:2511])
panGenome_nonORF <- panGenome_newRxn0[1:501,]
index1 <- which(panGenome_newRxn0$rxn %in% newOK ==TRUE)
newRxn_panGenome <- panGenome_nonORF[index1,]


#obtain the possible GPRs
unique_rxn <- unique(panGenome_newRxn0$formula)
index2 <- which(duplicated(panGenome_newRxn0$formula)==FALSE)
panGenome_newRxn_unique <- panGenome_newRxn0[index2,]
panGenome_newRxn_unique$query <- getMultipleReactionFormula(panGenome_newRxn0$query,panGenome_newRxn0$formula,panGenome_newRxn_unique$formula)







# panGenome annotation based on : http://eggnogdb.embl.de/#/app/home
egg1 <- scan("panGenome_annotation_egg/panProtein_s288_part1.txt.emapper.annotations", sep = "\n", what = "complex")
egg2 <- scan("panGenome_annotation_egg/panProtein_s288_part2.txt.emapper.annotations", sep = "\n", what = "complex")

egg10 <- read.delim2("panGenome_annotation_egg/panProtein_s288_part1.txt.emapper.annotations", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
egg20 <- read.delim2("panGenome_annotation_egg/panProtein_s288_part2.txt.emapper.annotations", header = FALSE, sep = "\t",stringsAsFactors = FALSE)

panGenome_eggnog <- rbind.data.frame(egg20, egg10)
colnames(panGenome_eggnog) <- c('query','seed_ortholog','evalue','score', 'predicted_name',
                             'GO_terms', 'kegg_ko', 'bigg_reactions', 'tax_scope', 'eggNOG_OGs',
                             'best_OG', 'COG_cat', 'eggNOG_HMM_desc')

## we only analysis the results mapping onto BiGG database
panGenome_eggnog0 <- select(panGenome_eggnog, query, bigg_reactions)
panGenome_eggnog0 <- filter(panGenome_eggnog0, bigg_reactions!="")
panGenome_eggnog0_rxn <- splitAndCombine(panGenome_eggnog0$bigg_reactions,panGenome_eggnog0$query,sep0 = ",")
colnames(panGenome_eggnog0_rxn) <- c('biggID','query')

#input bigg model reactions
bigg_models_reactions <- read_excel("BiGG/bigg_models_reactions.xlsx")

#find the reaction 
panGenome_eggnog0_rxn$name <- getMultipleReactionFormula(bigg_models_reactions$name,bigg_models_reactions$bigg_id,panGenome_eggnog0_rxn$biggID)
panGenome_eggnog0_rxn$reaction_string <- getMultipleReactionFormula(bigg_models_reactions$reaction_string,bigg_models_reactions$bigg_id,panGenome_eggnog0_rxn$biggID)
panGenome_eggnog0_rxn$database_links <- getMultipleReactionFormula(bigg_models_reactions$database_links,bigg_models_reactions$bigg_id,panGenome_eggnog0_rxn$biggID)


#find the new reaction
newBigg <- setdiff(panGenome_eggnog0_rxn$biggID[1:545], panGenome_eggnog0_rxn$biggID[545:1735])
panGenome_nonORF <- panGenome_newRxn0[1:501,]
index1 <- which(panGenome_newRxn0$rxn %in% newOK ==TRUE)
newRxn_panGenome <- panGenome_nonORF[index1,]
# we can't find the new reactions based on eggnog mapping onto the BiGG database











## analyse the eggnog mapping using the method_improved coverage but may lower the precision
# panGenome annotation based on : http://eggnogdb.embl.de/#/app/home
egg10 <- read.delim2("panGenome_annotation_egg/panProtein_s288_part1.txt.emapper.annotations_improved_coverage", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
egg20 <- read.delim2("panGenome_annotation_egg/panProtein_s288_part2.txt.emapper.annotations_improved_coverage", header = FALSE, sep = "\t",stringsAsFactors = FALSE)

panGenome_eggnog <- rbind.data.frame(egg20, egg10)
colnames(panGenome_eggnog) <- c('query','seed_ortholog','evalue','score', 'predicted_name',
                                'GO_terms', 'kegg_ko', 'bigg_reactions', 'tax_scope', 'eggNOG_OGs',
                                'best_OG', 'COG_cat', 'eggNOG_HMM_desc')

## we only analysis the results mapping onto BiGG database
panGenome_eggnog0 <- select(panGenome_eggnog, query, bigg_reactions)
panGenome_eggnog0 <- filter(panGenome_eggnog0, bigg_reactions!="")
panGenome_eggnog0_rxn <- splitAndCombine(panGenome_eggnog0$bigg_reactions,panGenome_eggnog0$query,sep0 = ",")
colnames(panGenome_eggnog0_rxn) <- c('biggID','query')

#input bigg model reactions
bigg_models_reactions <- read_excel("BiGG/bigg_models_reactions.xlsx")

#find the reaction 
panGenome_eggnog0_rxn$name <- getMultipleReactionFormula(bigg_models_reactions$name,bigg_models_reactions$bigg_id,panGenome_eggnog0_rxn$biggID)
panGenome_eggnog0_rxn$reaction_string <- getMultipleReactionFormula(bigg_models_reactions$reaction_string,bigg_models_reactions$bigg_id,panGenome_eggnog0_rxn$biggID)
panGenome_eggnog0_rxn$database_links <- getMultipleReactionFormula(bigg_models_reactions$database_links,bigg_models_reactions$bigg_id,panGenome_eggnog0_rxn$biggID)


#find the new reaction
newBigg <- setdiff(panGenome_eggnog0_rxn$biggID[1:545], panGenome_eggnog0_rxn$biggID[545:1735])

index1 <- which(panGenome_newRxn0$rxn %in% newOK ==TRUE)
newRxn_panGenome <- panGenome_nonORF[index1,]
# Still we can't find the new reactions based on eggnog mapping onto the BiGG database













































