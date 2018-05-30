library(readxl)
library(tidyverse)
library(stringr)
library(readr)
yeast_7_7  <- read_excel("yeast_7_7_corrected_gene_from_isce.xlsx")
###### function to splite all the related genes
getSpliteGene <- function(genelist, reactionName){
##split of or relation
#genelist <- yeast_7_7$Gene.reaction.association
#reactionName <-yeast_7_7$Rxn.name

rxnNum <- length(reactionName)


GR <- list() # gene relation
ss <- list()


Add_and <- function(x1){
  tt0 <- vector()
  tt <- x1
  tt_length <- length(tt)
  for (i in 1:tt_length){
    tt0[i] <- paste(i, tt[i], sep = ";")
  }
  return(tt0)
}


Add_or <- function(x1){ #x1 <- GR[[1]]
  tt0 <- vector()
  tt <- x1
  tt_length <- length(tt)
  for (i in 1:tt_length){
    tt0[i] <- paste(i, tt[i], sep = ";")
  }
  return(tt0)
}


## split the first order of gene relation

splitOrRelation <- function(genelist, reactionName){
  #genelist <- "((YAL023C and YDL095W) or YDL093W or YJR143C or YOR321W)"
  #reactionName <- "r_0016"
  GR <- list() # gene relation
  ss <- list()
if (!(str_detect(genelist, "\\) or \\(") | str_detect(genelist, "\\) or") | str_detect(genelist, "or \\(") | str_detect(genelist, "\\) and \\(") | str_detect(genelist, "\\) and") | str_detect(genelist, "and \\("))){
 
  GR[1] <- genelist
  ss[[1]] <- paste(reactionName,GR, sep = "@none;" )
  
} else {
  
  if (str_detect(genelist, "\\) or \\(")){
    GR[1] <- str_split(genelist,"\\) or \\(" )
    
  } else{
    GR[1] <- genelist
  }
  
  if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "\\) or")) {
    GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"\\) or")
  } else{
    GR[1] <- GR[1]
  }
  
  if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "or \\(")){
    GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"or \\(" )
   
  } else{
    GR[1] <- GR[1]
  }
  

  
  if (str_detect(genelist, "\\) and \\(")){
    GR[1] <- str_split(genelist,"\\) and \\(" )
    
  } else{
    GR[1] <- GR[1]
  }
  
  if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "\\) and")) {
    GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"\\) and")
  } else{
    GR[1] <- GR[1]
  }
  
  if (str_detect(paste0(unlist(GR[1]),collapse = "@@"), "and \\(")){
    GR[1] <- str_split(paste0(unlist(GR[1]),collapse = "@@"),"and \\(" )
    
  } else{
    GR[1] <- GR[1]
  }



##for genelist <- "((YAL023C and YDL095W) or YDL093W or YJR143C or YOR321W)" 
##reactionName <- "r_0016"
##tt <- GR[1]
  
  sl <-length(GR[[1]])
  for (i in 1:sl){
    if( (str_detect(GR[[1]][i],"\\(\\(")) | (str_detect(GR[[1]][i],"\\)\\)"))){
      GR[[1]][i] <- GR[[1]][i]
  } else if (str_detect(GR[[1]][i], "or")){
      GR[[1]][i]  <- str_split(GR[[1]][i], "or")
  } 
  }
  ##### above is for special condition
    
    
    
  GR[1] <- paste0(unlist(GR[1]),collapse = "@@")
  GR[1] <- str_split(GR[1],"@@")
  
  if (str_detect(genelist, "\\) or \\(") | str_detect(genelist, "\\) or") | str_detect(genelist, "or \\(")){
  ss[[1]] <- paste(reactionName,Add_or(GR[[1]]), sep = "@or" )
  } else {
  ss[[1]] <- paste(reactionName,Add_and(GR[[1]]), sep = "@and" )
  }
}

return(ss[[1]])

}

for (i in 1:rxnNum){
  
  ss[[i]] <- splitOrRelation(genelist[i],reactionName[i])

}


GR1 <- unlist(ss)

GR2 <-  str_replace_all(GR1,"\\(", "")
GR2 <-  str_replace_all(GR2,"\\)", "")
GR3 <- str_split(GR2,";")
tt0 <- length(GR3)
GR4 <- data.frame( ID= character(tt0),gene=character(tt0), stringsAsFactors = FALSE )
for (i in 1:tt0){
  GR4$ID[i] <- GR3[[i]][1]
  GR4$gene[i] <- GR3[[i]][2]
  
}

## split the first order of gene relation
rxnNum2 <- length(GR4$ID)
GR_2 <- list()
ss2 <- list()
for (i in 1:rxnNum2){
  if (str_detect(GR4$gene[i], " or ")){
    GR_2[i] <- str_split(GR4$gene[i]," or " )
    ss2[[i]] <- paste(GR4$ID[i],Add_or(GR_2[[i]]), sep = "@or" )
    
  } else if(str_detect(GR4$gene[i], " and ")){
    GR_2[i] <- str_split(GR4$gene[i]," and " )
    ss2[[i]] <- paste(GR4$ID[i],Add_and(GR_2[[i]]), sep = "@and" )
    
  } else{
    GR_2[i] <- GR4$gene[i]
    ss2[[i]] <- paste(GR4$ID[i],GR_2[[i]], sep = "@none;" )
    
  }
}


## obtain the final formula
GR_3 <- unlist(ss2)
GR_4 <- str_split(GR_3,";")
tt <- length(GR_4)
GR_5 <- data.frame( ID= character(tt),gene=character(tt), stringsAsFactors = FALSE )

for (i in 1:tt){
  GR_5$ID[i] <- GR_4[[i]][1]
  GR_5$gene[i] <- GR_4[[i]][2]
}

GR_or_and <- str_split(GR_5$ID,"@")

for (i in 1:tt){
  GR_5$IDnew[i] <- GR_or_and[[i]][1]
  GR_5$R1[i] <- GR_or_and[[i]][2]
  GR_5$R2[i] <- GR_or_and[[i]][3]
}

GPs_redesign <- select(GR_5,IDnew, R1, R2, gene) ##obtain the new format of GPRs, the can be base to add the new GPRs or correct GPRs

return(GPs_redesign)

}

yeast_7_7$Gene.reaction.association[is.na(yeast_7_7$Gene.reaction.association)] <- "NA"

GPs_redesign_yeast <- getSpliteGene(yeast_7_7$Gene.reaction.association, yeast_7_7$Rxn.name)


#process to combine all the splited GPR

# main function to combine all the splited GPR
getNewGPR <- function(GPs_redesign_yeast){
ND <- GPs_redesign_yeast

connect1 <- function(ss1_unique, gr2){
  if(str_detect(paste0(ss1_unique, collapse = " "),"and")){
    GR2 <- paste0(gr2,collapse = " and ")
    GR2 <- paste("(", GR2)
    GR2 <- paste(GR2, ")")
    
  } else if (str_detect(paste0(ss1_unique, collapse = " "),"or")){
    GR2 <- paste0(gr2,collapse = " or ") 
    GR2 <- paste("(", GR2)
    GR2 <- paste(GR2, ")")
    
  } else{
    GR2 <- gr2
  }
  
  return(GR2)
}

connect2 <- function(index2, ND){
if(str_detect(paste0(ND$R2[index2], collapse = " "),"and")){
    GR2 <- paste0(ND$gene[index2],collapse = " and ")
    GR2 <- paste("(", str_trim(GR2,"both"))
    GR2 <- paste(str_trim(GR2,"both"), ")")

} else if (str_detect(paste0(ND$R2[index2], collapse = " "),"or")){
    GR2 <- paste0(ND$gene[index2],collapse = " or ") 
    GR2 <- paste("(", str_trim(GR2,"both"))
    GR2 <- paste(str_trim(GR2,"both"), ")")
    
} else{
    GR2 <- ND$gene[index2]
}

  return(GR2)
}

getGPR <- function(rxn_list, one_unique_rxn){
  index1 <- which(rxn_list %in% one_unique_rxn ==TRUE)
  ss1 <- rxn_sort[index1]
  ss1_unique <- unique(ss1)
  tt <-length(ss1_unique)
  gr2 <- vector()
  index2 <- list()
  for (i in 1:tt) {
    index2[[i]] <- which(rxn_sort %in% ss1_unique[i])
    gr2[i] <- connect2(index2[[i]],ND)
  }
  
  GPR <- connect1(ss1_unique, gr2)
}



rxn <-ND$IDnew
rxn_unique <- unique(rxn)
rxn_sort <- paste0(ND$IDnew,ND$R1)
tt <- length(rxn_unique)
newGPR <- data.frame(ID=character(tt),stringsAsFactors = FALSE)
newGPR$ID <- rxn_unique
for (j in 1:tt){
  newGPR$GPR[j] <- getGPR(rxn,rxn_unique[j])
}
newGPR$GPR <- str_replace_all(newGPR$GPR, "\\( ", "\\(")
newGPR$GPR <- str_replace_all(newGPR$GPR, " \\)", "\\)")

return(newGPR)

}

# some preprocess before combine: remove the space
GPs_redesign_yeast$IDnew <- str_trim(GPs_redesign_yeast$IDnew,"both")
GPs_redesign_yeast$R1 <- str_trim(GPs_redesign_yeast$R1,"both")
GPs_redesign_yeast$R2 <- str_trim(GPs_redesign_yeast$R2,"both")
GPs_redesign_yeast$gene <- str_trim(GPs_redesign_yeast$gene,"both")

# use function to combine all the splited GPR
newGPR <- getNewGPR(GPs_redesign_yeast)


##test the resuls
yeast_7_7$NEWGPR <- newGPR$GPR
yeast_7_7$comparison <- newGPR$GPR== yeast_7_7$Gene.reaction.association
wrong <- filter(yeast_7_7, comparison==FALSE)
wrong$NEWGPR
wrong$Gene.reaction.association
wrong$Rxn.name
