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
panGenome_kegg <- read.delim2("panGenome_annotation from eggnog and kegg/panProtein_annotation_kegg.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(panGenome_kegg) <- c('query','ko')
panGenome_kegg0 <- filter(panGenome_kegg, ko != "")


#input the KO annotation from kegg database
#KO_rxn mapping
KO_rxn <- read.delim2("panGenome_annotation from eggnog and kegg/reaction-KO from kegg .txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(KO_rxn) <- c('rxnID','ko')
KO_rxn$ko <- str_replace_all(KO_rxn$ko, "ko:", "")
# rxn formual in kegg
rxn_kegg <- read.delim2("panGenome_annotation from eggnog and kegg/reaction_kegg summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#KO_pathway mapping
KO_pathway <- read.delim2("panGenome_annotation from eggnog and kegg/ko_pathway.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(KO_pathway) <- c('pathway','ko')
KO_pathway$ko <- str_replace_all(KO_pathway$ko, "ko:", "")
KO_pathway <- filter(KO_pathway,str_detect(KO_pathway$pathway,"ko")==FALSE)
#pathway annotation
pathway <- read.delim2("panGenome_annotation from eggnog and kegg/pathway_list_kegg.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(pathway) <- c('pathway','pathway_name')


#get the rxn based on KO
panGenome_kegg0$rxns <- getMultipleReactionFormula(KO_rxn$rxnID,KO_rxn$ko,panGenome_kegg0$ko)
panGenome_newRxn <- filter(panGenome_kegg0,!is.na(rxns))

panGenome_kegg <- splitAndCombine(panGenome_newRxn$rxns, panGenome_newRxn$query, sep0 = ";")
colnames(panGenome_kegg) <- c('rxn','query')

panGenome_kegg$formula <- getMultipleReactionFormula(rxn_kegg$reaction,rxn_kegg$rxnID,panGenome_kegg$rxn)
panGenome_kegg$description <- getMultipleReactionFormula(rxn_kegg$name,rxn_kegg$rxnID,panGenome_kegg$rxn)
panGenome_kegg$ko <- getMultipleReactionFormula(panGenome_kegg0$ko,panGenome_kegg0$query,panGenome_kegg$query)
panGenome_kegg$pathway <- getMultipleReactionFormula(KO_pathway$pathway,KO_pathway$ko,panGenome_kegg$ko)

gene_pathway <- splitAndCombine(panGenome_kegg$pathway, panGenome_kegg$query, sep0 = ";")
colnames(gene_pathway)  <- c('pathID','gene')
gene_pathway$pathway_name <- getMultipleReactionFormula(pathway$pathway_name,pathway$pathway,gene_pathway$pathID)
#change the pathway into sce
gene_pathway$pathID <- str_replace_all(gene_pathway$pathID,"path:map","sce")
gene_pathway$pathway_name0 <- paste(gene_pathway$pathID, gene_pathway$pathway_name, sep = "  ")
#merge the pathways
panGenome_kegg$pathway_name0 <- getMultipleReactionFormula(gene_pathway$pathway_name0,gene_pathway$gene,panGenome_kegg$query)

#obtain the possible GPRs
unique_rxn <- unique(panGenome_kegg$formula)
index2 <- which(duplicated(panGenome_kegg$formula)==FALSE)
panGenome_Rxn_kegg <- panGenome_kegg[index2,]
panGenome_Rxn_kegg$query <- getMultipleReactionFormula(panGenome_kegg$query,panGenome_kegg$formula,panGenome_Rxn_kegg$formula)
write.table(panGenome_Rxn_kegg, 'panGenome_Rxn_kegg.txt', row.names = FALSE, sep = '\t')

#standization of reactions using Rhea and metnetx database
#find the rheaID and chebiID based on rxnID
#find the chebiID based on rheaID
#find the stardard metabolite forumula based on chebiID or rheaID
#this step will use a standard pipeline


#input the standard formula of reacions and metabolites
#obtain the new Rxns with the standard formula of reacions and metabolites
newRxnID <- setdiff(panGenome_kegg$rxn[1:501], panGenome_kegg$rxn[502:2511])
index1 <- which(panGenome_Rxn_kegg$rxn %in% newRxnID ==TRUE)
newRxn_panGenome <- panGenome_Rxn_kegg[index1,]






# panGenome annotation based on : http://eggnogdb.embl.de/#/app/home
# it is quite detailed for the panGenome annotation from eggnog web service
# the results could the direct evidences for the article
# with the mapping onto BiGG, no new reactions were found for these new panID compared with mapping onto kegg database
# it may be due to the fact that the metabolic coverage in the bigg database is quite small
egg1 <- scan("panGenome_annotation from eggnog and kegg/panProtein_s288_part1.txt.emapper.annotations", sep = "\n", what = "complex")
egg2 <- scan("panGenome_annotation from eggnog and kegg/panProtein_s288_part2.txt.emapper.annotations", sep = "\n", what = "complex")

egg10 <- read.delim2("panGenome_annotation from eggnog and kegg/panProtein_s288_part1.txt.emapper.annotations", header = FALSE, sep = "\t",stringsAsFactors = FALSE)
egg20 <- read.delim2("panGenome_annotation from eggnog and kegg/panProtein_s288_part2.txt.emapper.annotations", header = FALSE, sep = "\t",stringsAsFactors = FALSE)

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
panGenome_eggnog0_rxn$biggID <- str_trim(panGenome_eggnog0_rxn$biggID, side = "both")
newBigg <- setdiff(panGenome_eggnog0_rxn$biggID[1:545], panGenome_eggnog0_rxn$biggID[545:1735])

#however we can't find the new reactions based on eggnog mapping onto the BiGG database












































