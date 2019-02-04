# input the main_function and package
source('main_function_map.R')


########small task###update the metabolites formula
#read the model

rxn <- read_excel("data/yeastGEM_latest version.xls", 
                                             sheet = "Reaction List")
metabolite <-  read_excel("data/yeastGEM_latest version.xls", 
                          sheet = "Metabolite List")
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

rxn_split <- splitAndCombine0(rxn$Reaction,rxn$Abbreviation,sep=" ")
rxn_split$v3 <- getSingleReactionFormula(metabolite$`Metabolite description`,metabolite$`Metabolite name`,rxn_split$v1)

for (i in 1:length(rxn_split$v2)){
  if(rxn_split$v3[i]=="NA"){
    rxn_split$v3[i] <- rxn_split$v1[i]
  } else{
    rxn_split$v3[i] <- rxn_split$v3[i]
  }
  

}
rxn$Description_new <- getMultipleReactionFormula(rxn_split$v3,rxn_split$v2,rxn$Abbreviation)
rxn$Description_new <- str_replace_all(rxn$Description_new,";"," ")
write.table(rxn, "result/yeastGEM with update reaction formula.txt", row.names = FALSE, sep = "\t")
########small task###update the metabolites formula

########small task 2###establish the relation between rxn
 # obtain the rxnID and metabolites
index <- which(str_detect(rxn_split$v3,"\\[")==TRUE)
rxn_metabolite <- rxn_split[index,]

 # find the currency metabolites
rxn_metabolite_withoutCompartment <- rxn_metabolite
rxn_metabolite_withoutCompartment$v3 <- str_replace_all(rxn_metabolite$v3,"\\[.*?\\]","")
analysis_metabolites <- rxn_metabolite_withoutCompartment %>% 
  count(v3) %>% ## calculate the number of each metabolite
  arrange(., desc(n)) ## order the metabolites based on the number
currency_metabolites <- analysis_metabolites$v3[1:14] ## find the currency metabolites


  # remove the currency metabolites 
index0 <- which(rxn_metabolite_withoutCompartment$v3 %in% currency_metabolites==FALSE)
rxn_metabolite_new <- rxn_metabolite[index0,]
colnames(rxn_metabolite_new) <- c("metabolitesID","rxn","metabolites")



getConnectedReaction <- function(ss0,m,rxn_metabolite_new,rxn_list){ 
  ## ss0: a rxn id list
  ## m: represent the position of aimed rxn id, when establish the relation
index1 <- which(rxn_metabolite_new$rxn %in% ss0==TRUE)
m0 <- rxn_metabolite_new$metabolites[index1]
rxn_list <- unique(rxn_metabolite_new$rxn)

rx_list_remaining <- rxn_list[-(1:m)] ## to remove the reactions before the m'th rxn
index2 <- which(rxn_metabolite_new$rxn %in% rx_list_remaining==TRUE)
rxn_metabolite_remaining <- rxn_metabolite_new[index2,]

index3 <- which(rxn_metabolite_remaining$metabolites %in% m0 ==TRUE)
rxn_connected <- unique(rxn_metabolite_remaining$rxn[index3])
rxn_connected <- paste(rxn_connected,collapse = ";")

return(rxn_connected)

}

# test example of this function:
# getConnectedReaction(rxn_list[3],m=3,rxn_metabolite_new,rxn_list)


# define a data.fram to establish the reactions relation
rxn_list <- unique(rxn_metabolite_new$rxn) ## get all the unique reaction ID
reaction_reaction <- data.frame( rxnID = character(length(rxn_list)), stringsAsFactors = FALSE)
reaction_reaction$rxnID <- rxn_list

for (i in 1:(length(rxn_list)-1)){
  
  reaction_reaction$connected_rxn[i] <- getConnectedReaction(rxn_list[i],m=i,rxn_metabolite_new,rxn_list)
  
}


index4 <- which (reaction_reaction$connected_rxn != "")
reaction_reaction_final <- reaction_reaction[index4,]
reaction_reaction_final0 <- splitAndCombine0(reaction_reaction_final$connected_rxn,reaction_reaction_final$rxnID,sep=";")
reaction_reaction_final0 <- select(reaction_reaction_final0, v2,v1)
colnames(reaction_reaction_final0) <- c('from','to')
########small task 2###establish the relation between rxn



#classify the reactions based on subsystems
rxn_system <- select(rxn,Abbreviation,Subsystem_new)
index00 <- which(str_detect(rxn_system$Subsystem_new,"transport")==TRUE)
rxn_system$Subsystem_new[index00] <-"transport"
subsystem <- rxn_system %>% 
  count(Subsystem_new) %>%
  arrange(.,desc(n))

# display the reaction relation from certain subsystems
choosed_subsystem <- c(subsystem$Subsystem_new[3])
index_s1 <- which(rxn_system$Subsystem_new %in% choosed_subsystem ==TRUE)
rxn_s1 <- rxn_system$Abbreviation[index_s1]

s1 <- which(reaction_reaction_final0$from %in% rxn_s1 ==TRUE)
s2 <- which(reaction_reaction_final0$to %in% rxn_s1 ==TRUE)
s3 <- intersect(s1,s2)
links_s1 <- reaction_reaction_final0[s3,]

nodes0 <- rxn_s1
rxnLinks <- links_s1
colnames(rxnLinks) <- c("source0","target0")
rxnLinks$value <- 0.1

rxnNodes <- data.frame(name=character(length(nodes0)),group=character(length(nodes0)),size=character(length(nodes0)),stringsAsFactors = FALSE)
rxnNodes$name <- nodes0
rxnNodes$size <- 8
rxnNodes$order <- (0:(length(nodes0)-1))

#give subsystem for the reactions
rxnNodes$group <- getSingleReactionFormula(rxn_system$Subsystem_new,rxn_system$Abbreviation,rxnNodes$name)
rxnLinks$source <- as.numeric(getSingleReactionFormula(rxnNodes$order,rxnNodes$name,rxnLinks$source0))
rxnLinks$target <- as.numeric(getSingleReactionFormula(rxnNodes$order,rxnNodes$name,rxnLinks$target0))

rxnLinks0 <- select(rxnLinks, source, target, value)
rxnNodes0 <- select(rxnNodes,name, group, size)

# simpleNetwork function only show the connected nodes while forceNetwork could show all the nodes
simpleNetwork(links_s1, fontSize = 18, opacity = 1, zoom = TRUE,fontFamily = "sans-serif")
# forceNetwork more detailed
forceNetwork(Links = rxnLinks0, Nodes = rxnNodes0, Source = "source",
             Target = "target", Value = "value", NodeID = "name",Nodesize = "size", radiusCalculation = "Math.sqrt(d.nodesize)", fontSize = 12,
             Group = "group", opacity = 1, zoom = TRUE, bounded = FALSE,
             fontFamily = "serif", opacityNoHover = 1, legend=T)





### display the whole graph
links0 <- reaction_reaction_final0[1:32191,]
nodes0 <- unique(c(links0$from, links0$to))
#using the data for the reactions data
rxnLinks <- links0
colnames(rxnLinks) <- c("source0","target0")
rxnLinks$value <- 0.1

rxnNodes <- data.frame(name=character(length(nodes0)),group=character(length(nodes0)),size=character(length(nodes0)),stringsAsFactors = FALSE)
rxnNodes$name <- nodes0
rxnNodes$group <- 1
rxnNodes$size <-8
rxnNodes$order <- (0:(length(nodes0)-1))


rxnLinks$source <- getSingleReactionFormula(rxnNodes$order,rxnNodes$name,rxnLinks$source0)
rxnLinks$target <- getSingleReactionFormula(rxnNodes$order,rxnNodes$name,rxnLinks$target0)

rxnLinks0 <- select(rxnLinks, source, target, value)
rxnNodes0 <- select(rxnNodes,name, group, size)
forceNetwork(Links = rxnLinks0, Nodes = rxnNodes0, Source = "source",
             Target = "target", Value = "value", NodeID = "name",Nodesize = "size", radiusCalculation = "Math.sqrt(d.nodesize)",fontSize = 12,
             Group = "group", opacity = 1, zoom = TRUE, bounded = FALSE,
             fontFamily = "serif", opacityNoHover = 1)
simpleNetwork(links0, fontSize = 18, opacity = 1, zoom = TRUE,fontFamily = "sans-serif")







