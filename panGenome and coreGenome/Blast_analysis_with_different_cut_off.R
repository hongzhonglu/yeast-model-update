library(stringr)
library(tidyverse)
library(ggplot2)
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

#part1
#input the data blast with uniprot
u_rxn_45_100 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_45_100.txt', stringsAsFactors = FALSE)
u_rxn_35_100 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_35_100.txt', stringsAsFactors = FALSE)
u_rxn_55_100 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_55_100.txt', stringsAsFactors = FALSE)
u_rxn_65_100 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_65_100.txt', stringsAsFactors = FALSE)
u_rxn_75_100 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_75_100.txt', stringsAsFactors = FALSE)
#u_rxn_85_100 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_85_100.txt', stringsAsFactors = FALSE)
u_rxn_45_50 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_45_50.txt', stringsAsFactors = FALSE)
u_rxn_45_150 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_45_150.txt', stringsAsFactors = FALSE)
u_rxn_45_200 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_45_200.txt', stringsAsFactors = FALSE)
u_rxn_45_250 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_45_250.txt', stringsAsFactors = FALSE)

dfList <- list(u_rxn_45_100,u_rxn_35_100,u_rxn_55_100,u_rxn_65_100,u_rxn_75_100,u_rxn_45_50,u_rxn_45_150,u_rxn_45_200,u_rxn_45_250)

calUniqueNum =function (ss){
  #test
  #ss <- u_rxn_35_100$panID
  ss0 <- str_split(ss, ";")
  ss1 <- unlist(ss0)
  tt <- length(unique(ss1))
  return(tt)
}


condition_list <- c('u_rxn_45_100','u_rxn_35_100','u_rxn_55_100','u_rxn_65_100','u_rxn_75_100','u_rxn_45_50','u_rxn_45_150','u_rxn_45_200','u_rxn_45_250')
rxn_sum <- data.frame(condition =condition_list, stringsAsFactors = FALSE)

for (i in seq(length(dfList))){
  ss <- dfList[[i]]
  rxn_sum$rxn_num[i] <- calUniqueNum(ss$panID)
    
}


#effects of idendity
rxn_sum_fixed_bitscore <- rxn_sum[str_detect(rxn_sum$condition,"_100"),]
rxn_sum_fixed_bitscore$condition <- c('45%','35%','55%','65%','75%')

# Basic barplot
ggplot(data=rxn_sum_fixed_bitscore, aes(x=condition, y=rxn_num)) +
  geom_bar(stat="identity") +
  xlab("Identity") + ylab("Count of new proteins")



#effects of bitscore
rxn_sum_fixed_identity <- rxn_sum[str_detect(rxn_sum$condition,"_45"),]
rxn_sum_fixed_identity$condition <- c('100','50','150','200','250')

ggplot(data=rxn_sum_fixed_identity, aes(x = reorder(condition, -rxn_num), y = rxn_num)) + 
  geom_bar(stat = "identity") +
  xlab("Bitscore") + ylab("Count of new proteins")



#part2
#input the data blast with bigg
b_rxn_45_100 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_45_100.txt', stringsAsFactors = FALSE)
b_rxn_55_100 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_55_100.txt', stringsAsFactors = FALSE)
b_rxn_65_100 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_65_100.txt', stringsAsFactors = FALSE)
b_rxn_75_100 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_75_100.txt', stringsAsFactors = FALSE)
#u_rxn_85_100 <- read.delim2('blast with uniprot using different cut-off/new rxn from pan based on uniprot_85_100.txt', stringsAsFactors = FALSE)
b_rxn_35_100 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_35_100.txt', stringsAsFactors = FALSE)
b_rxn_45_150 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_45_150.txt', stringsAsFactors = FALSE)
b_rxn_45_200 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_45_200.txt', stringsAsFactors = FALSE)
b_rxn_45_250 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_45_250.txt', stringsAsFactors = FALSE)
b_rxn_45_50 <- read.delim2('blast with bigg using different cut-off/new rxn from pan_45_50.txt', stringsAsFactors = FALSE)


dfList <- list(b_rxn_45_100,b_rxn_55_100,b_rxn_65_100,b_rxn_75_100,b_rxn_35_100,b_rxn_45_150,b_rxn_45_200,b_rxn_45_250,b_rxn_45_50)


condition_list <- c('b_rxn_45_100','b_rxn_55_100','b_rxn_65_100','b_rxn_75_100','b_rxn_35_100','b_rxn_45_150','b_rxn_45_200','b_rxn_45_250','b_rxn_45_50')
rxn_sum <- data.frame(condition =condition_list, stringsAsFactors = FALSE)

for (i in seq(length(dfList))){
  ss <- dfList[[i]]
  rxn_sum$rxn_num[i] <- calUniqueNum(ss$panID)
  
}


#effects of idendity
rxn_sum_fixed_bitscore <- rxn_sum[str_detect(rxn_sum$condition,"_100"),]
rxn_sum_fixed_bitscore$condition <- c('45%','55%','65%','75%','35%')

# Basic barplot
ggplot(data=rxn_sum_fixed_bitscore, aes(x=condition, y=rxn_num)) +
  geom_bar(stat="identity") +
  xlab("Identity") + ylab("Count of new proteins")



#effects of bitscore
rxn_sum_fixed_identity <- rxn_sum[str_detect(rxn_sum$condition,"_45"),]
rxn_sum_fixed_identity$condition <- c('100','150','200','250','50')

ggplot(data=rxn_sum_fixed_identity, aes(x = reorder(condition, -rxn_num), y = rxn_num)) + 
  geom_bar(stat = "identity") +
  xlab("Bitscore") + ylab("Count of new proteins")

