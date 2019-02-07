library(stringr)
library(tidyverse)
library(readxl)
library(igraph)
library(networkD3)
library(hongR)


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


rxn_metabolite_mapping <- select(rxn,Abbreviation, Reaction)
rxn_metabolite_mapping$Reaction <- str_replace_all(rxn_metabolite_mapping$Reaction, "\\-\\>","\\<\\=\\>" )


getMappingReaction <- function(reaction,id) {
  # input: reaction equation
  # input: reaction id
  # this function is used to establish the relation between metabolite according to certain reaction
  
rxn_metabolite_split0 <- splitAndCombine0(reaction,id,sep=" \\<\\=\\> ")
reactants <- str_split(rxn_metabolite_split0$v1[1], " ")
reactants <- unlist(reactants)
reactants_index <- which(str_detect(reactants,"\\[")==TRUE)
reactants <- reactants[reactants_index]

products <- str_split(rxn_metabolite_split0$v1[2]," ")
products <- unlist(products)
products_index <- which(str_detect(products,"\\[")==TRUE)
products <- products[products_index]


m <- length(reactants)
n <- length(products)
mn <- c()
for (i in seq(m)){
  for(j in seq(n)){
    mn <- c(paste(reactants[i],products[j],sep = ";"),mn)
  }
}

return(mn)
}


## process to establish the metabolite mapping relation
ss <- list()
ss0 <- list()
for (i in 1:3496){
ss[[i]] <- getMappingReaction(rxn_metabolite_mapping$Reaction[i],rxn_metabolite_mapping$Abbreviation[i])
ss0[[i]] <- paste(rxn_metabolite_mapping$Abbreviation[i],ss[[i]],sep = ";" )
}

mappedMetabolites <- unlist(ss0)
head(mappedMetabolites)

ss1 <- str_split(mappedMetabolites,";")
ss2 <- data.frame(id = character(length(ss1)), stringsAsFactors = FALSE)
for (i in seq(length(ss1))){
  ss2$id[i] <- ss1[[i]][1]
  ss2$from[i] <- ss1[[i]][2]
  ss2$to[i] <- ss1[[i]][3]
  ss2$combine[i] <- paste(ss1[[i]][2],ss1[[i]][3])
}

ss2 <- ss2[which(ss2$to!=""),]

ss3 <- ss2[!duplicated(ss2$combine),] # remove the duplicated relation of metabolites ?

ss3 <- select(ss3, from, to)

simpleNetwork(ss3[1:100,], fontSize = 18, opacity = 1, zoom = TRUE,fontFamily = "sans-serif")

# remove the currency metabolite
ss4_1 <- select(ss3,from)
colnames(ss4_1) <- "metabolite"
ss4_2 <- select(ss3, to)
colnames(ss4_2) <- "metabolite"
ss4 <-rbind(ss4_1,ss4_2)
ss4$name <- getSingleReactionFormula(metabolite$`Metabolite description`,metabolite$`Metabolite name`,ss4$metabolite)
ss4$name_simple <-  str_replace_all(ss4$name,"\\[.*?\\]","")

metabolite_analysis <- ss4 %>%
  count(name_simple) %>%
  arrange(.,desc(n))

currency_metabolite <- metabolite_analysis$name_simple[1:14]
index1 <- which(ss4$name_simple %in% currency_metabolite == TRUE)

# find the related id of currency metabolite
id <- unique(ss4$metabolite[index1])

# remove the currency metabolite from links dataframe
index2_1 <- which(ss3$from %in% id ==FALSE)
index2_2 <- which(ss3$to %in% id ==FALSE)
index2 <- intersect(index2_1, index2_2)

ss30 <- ss3[index2,]



require(igraph)
library(tidyverse)
library(stringr)

links <- ss30
links$weight <- 10
nodes0 <- unique(c(ss30$from,ss30$to))
nodes <- data.frame(id=nodes0,stringsAsFactors = FALSE)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 



#Breadth-first search
#s_0563[c] D-glucose[C]
ss <- bfs(net, root= "s_0563[c]", "out",unreachable = FALSE,
          order=TRUE, rank=TRUE, father=TRUE, pred=TRUE,
          succ=TRUE, dist=TRUE)


node_parameter <- select(nodes, id)
node_parameter$dist <- ss$dist

s1 <- which(node_parameter$dist %in% 1 ==TRUE)
node_parameter[s1,]


#It is very slow using the followed two codes
#all_simple_paths(net,from="s_0563[c]", to="s_1405[c]")
#shortest_path <- all_shortest_paths(net,from="s_0563[c]", to="s_1405[c]")



# this two functions could find all the children node based on input node
getNextLayer <- function(NET=net, NODE = nodes, root ){

  ss <- bfs(NET, root, "out", unreachable = FALSE,
            order=TRUE, rank=TRUE, father=TRUE, pred=TRUE,
            succ=TRUE, dist=TRUE)
  
  node_parameter <- select(NODE, id)
  node_parameter$dist <- ss$dist
  node1 <- filter(node_parameter, dist==1)
  node2 <- node1$id
  return(node2)
}

getMultiNodes <- function (multiRoot){
  #mutiRoot = tt[[5]]
  ss <- unlist(multiRoot)
  nn <- list()
  mm <- length(ss)
  for (i in 1:mm){
    nn[[i]] <- getNextLayer(NET = net, NODE = nodes, ss[i] )
    
  }
  
  tt <- unique(unlist(nn))
  return(tt)
}



#function get the connection

getPath <- function (t){
  tt2 <- t
  nn <- length(tt2)
  length0 <- vector()
  n0 <- list()
  for (i in 1:nn){
    n0[[i]] <- getNextLayer(root = tt2[[1]][length(tt2[[1]])])
    length0[i] <- length(n0[[i]])
  }
  
  #next layer
  nn0 <- unlist(n0)
  
  tt20 <- list()
  for (i in 1:nn){
    tt20[[i]] <- paste0(tt2[[i]], collapse = "@")
  }
  
  mm0 <- list()
  
  for (i in 1:nn){
    mm0[[i]] <- replicate(length0[i], tt20[[i]])
  }
  
  mm1 <- unlist(mm0)
  
  # combine the current layer and next layer
  total <- sum(length0)
  
  mm1 <- paste(mm1,nn0, sep = "@")  
  
  
  mm2 <- str_split(mm1,"@")
  
  return(mm2)
  
}




tt1 <- list()
tt1[[1]] <- "s_0563[c]"
tt2 <- getPath(tt1)
tt3 <- getPath(tt2)
tt4 <- getPath(tt3)
tt5 <- getPath(tt4)
tt6 <- getPath(tt5)
tt7 <- getPath(tt6)

# It is become very slow after the 7th step. For more step, new method will be used
#tt8 <- getPath(tt7)
#tt9 <- getPath(tt8)



## function get connection using unique metabolite to accelerate this process
getUniqueMetabolitePath <- function(t){
tt4 <- t
nn4 <- length(tt4)
mm4 <- vector()

for (i in 1:nn4){
mm4[i] <- tt4[[i]][length(tt4[[i]])]
}

mm40 <- unique(mm4)
tt40 <- list()
for (i in 1:length(mm40)){
  
  tt40[[i]] <- mm40[i]
  
}
tt50 <- getPath(tt40)

return(tt50)
}

tt8 <- getUniqueMetabolitePath(tt7)
tt9 <- getUniqueMetabolitePath(tt8)
tt10 <- getUniqueMetabolitePath(tt9)
tt11 <- getUniqueMetabolitePath(tt10)
tt12 <- getUniqueMetabolitePath(tt11)
tt13 <- getUniqueMetabolitePath(tt12)
tt14 <- getUniqueMetabolitePath(tt13)
tt15 <- getUniqueMetabolitePath(tt14)




##example

#s_1405[c]	riboflavin[c]

met <- "s_1405[c]"
met %in% unique(unlist(tt3))

# analysis
ll <- length(tt3)
hh <- vector()
for (i in 1:ll){
  hh[i] <- tt3[[i]][3]
  
}

which(hh %in% met ==TRUE)

tt3[7322]

index0 <- which(ss30$from %in% "s_0027[m]" ==TRUE)

ss30[index0,]
