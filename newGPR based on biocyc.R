library(readxl)
library(stringr)
library(tidyverse)
reaction <- read_excel("reaction.xlsx")
reaction$X__1 <- str_replace_all(reaction$X__1, "↔", "<=>")
reaction$X__1 <- str_replace_all(reaction$X__1, " → ", " => ")  ## some metabolites have "→"
reaction$X__1 <- str_replace_all(reaction$X__1, "→", "->")  ## some metabolites have "→"
reaction$X__1 <- str_replace_all(reaction$X__1, "α", "alpha")
reaction$X__1 <- str_replace_all(reaction$X__1, "β", "beta")
reaction$X__1 <- str_replace_all(reaction$X__1, "ω", "omega")

#function 1
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

#function 2 estimate whether the gene exist in presnet model

geneExist <- function (original, newgenelist){  ##input: original is from yeast model v7.7; 
  index <- vector()
  for (i in 1: length(newgenelist)){
    if (length(which(newgenelist[i] %in% original))){
      index[i] <- "YES"
    } else {
      index[i] <- "NO"
    }
  }
  return(index)
}


## Obtain the sytematic name of gene related with the reactions
gene_standard_name <- read_excel("gene_name.xlsx")
gene_standard_name <- select(gene_standard_name, `Accession-1`, `Common-Name`)
colnames(gene_standard_name)<- c("systematic_name","comman_name")
genelist_v7_7 <- read.csv("genelist in v7.7.csv", sep = ";", stringsAsFactors = FALSE)


## Establish relations between genes and reactions one by one
index_r <- which(!is.na(reaction$X__1))

reaction1 <- reaction
for (i in 2:3036) {
  if(is.na(reaction1$X__1[i])){
    reaction1$X__1[i] <- reaction1$X__1[i-1]
  } else{
    reaction1$X__1[i] <- reaction1$X__1[i]
  }
  
} ## in reaction1, for each gene, a reaction can be found

rxn_withGENE <- select(reaction1, X__1, gene) 
rxn_withGENE <- filter(reaction1, !is.na(gene))
GR <- select(rxn_withGENE, gene, X__1)
GR$sytematic_name <- getMultipleReactionFormula(gene_standard_name$systematic_name,gene_standard_name$comman_name,GR$gene)
GR0 <- select(GR, X__1, sytematic_name)

splitAndCombine <- function(gene, rxn) { ##one rxn has several genes, this function was used to splite the genes

   gene <- str_split(gene,";")
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

GR_biocyc <- splitAndCombine(GR0$sytematic_name, GR0$X__1)
colnames(GR_biocyc) <- c("gene","reaction_biocyc")
GR_biocyc$sign <- geneExist(genelist_v7_7$geneNames, GR_biocyc$gene)
newGR_biocyc <- filter(GR_biocyc,sign =="NO")
write.table(newGR_biocyc, "newGR_biocyc.txt", row.names = FALSE, sep = "\t")
###########################################################################################################newGPR obtained from biocyc

## gene annotation ec
## gene annotation from biocyc
gene_ec_biocyc <- read.table("gene_ec_biocyc.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  select(., Gene.Name, Reactions.of.gene) %>%
  filter(., str_detect(Reactions.of.gene, "\\."))
gene_ec_biocyc$Reactions.of.gene <- str_replace_all(gene_ec_biocyc$Reactions.of.gene,"-RXN","")
gene_ec_biocyc$Reactions.of.gene <- str_replace_all(gene_ec_biocyc$Reactions.of.gene,"RXN0-","")
gene_ec_biocyc$Reactions.of.gene <- str_replace_all(gene_ec_biocyc$Reactions.of.gene,"RXN-","")

splitAndCombine <- function(gene, rxn) { ##one rxn has several genes, this function was used to splite the genes
  
  gene <- str_split(gene,"//")
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

gene_ec_biocyc0 <- splitAndCombine(gene_ec_biocyc$Reactions.of.gene, gene_ec_biocyc$Gene.Name) %>%
  filter(., str_detect(v1,"\\."))


## get annotation and subsystem for each gene
gene_pathway_biocyc <- read.table("All_genes_pathway_biocyc.txt", header = TRUE, sep="\t", stringsAsFactors = FALSE)
gene_pathway_biocyc$gene_standard_name <- getMultipleReactionFormula(gene_standard_name$systematic_name,gene_standard_name$comman_name,gene_pathway_biocyc$Gene.Name)
gene_pathway_biocyc$check <- gene_pathway_biocyc$gene_standard_name == gene_pathway_biocyc$Accession.1
write.table(gene_pathway_biocyc, "gene_pathway and annotation in biocyc.txt", row.names = FALSE, sep="\t") ## >>> for manual check
##>>>>>after manual check
gene_pathway_biocyc0 <- read_excel("gene_pathway and annotation in biocyc after manual check.xlsx", 
                      sheet = "Sheet1")

gene_pathway_biocyc0$ec <- getMultipleReactionFormula(gene_ec_biocyc0$v1,gene_ec_biocyc0$v2,gene_pathway_biocyc0$Gene.Name)
gene_pathway_biocyc0$reaction_biocyc <- getMultipleReactionFormula(GR_biocyc$reaction_biocyc,GR_biocyc$gene,gene_pathway_biocyc0$gene_standard_name)
write.table(gene_pathway_biocyc0, "gene_pathway and annotation in biocyc0.txt", row.names = FALSE, sep="\t")


##get reactions according to EC
## get reactions from Rhea
#reaction_rhea<- read.csv("rhea reaction summary.csv", sep = ",", stringsAsFactors = FALSE)
#newGP$reaction_R <- getMultipleReactionFormula(reaction_rhea$formula,reaction_rhea$EC,newGP$ec)
#newGP$keggID_R <- getMultipleReactionFormula(reaction_rhea$keggID,reaction_rhea$EC,newGP$ec)
#newGP$masterID_R <- getMultipleReactionFormula(reaction_rhea$masterID,reaction_rhea$EC,newGP$ec)

##get reaction from Brenda
#reaction_brenda <- read.csv("Reactions_BKMS.csv", sep = ";", stringsAsFactors = FALSE)
#newGP$reaction_B <- getMultipleReactionFormula(reaction_brenda$Reaction,reaction_brenda$EC.Number,newGP$ec)
#newGP$keggID_B <- getMultipleReactionFormula(reaction_brenda$Reaction.ID.KEGG,reaction_brenda$EC.Number,newGP$ec)
#newGP$keggPathwayName_B <- getMultipleReactionFormula(reaction_brenda$KEGG.Pathway.Name,reaction_brenda$EC.Number,newGP$ec)
#newGP$brendaID_B <- getMultipleReactionFormula(reaction_brenda$Reaction.ID.BRENDA,reaction_brenda$EC.Number,newGP$ec)
#newGP0 <- filter(newGP, reaction_R !="NA" | reaction_B != "NA")


#write.table(newGP0, "newGP_biocyc.txt", row.names = FALSE, sep = "\t")






