##txt6 is the main reaction format information from xml of cellDesigner
##txt6 <- readLines(file("template6")) # repeated part: information of each reactions
##split txt6_repeated_part into several part
#R1 <- readLines(file("txt6_repeated_part1")) # define the reaction metaid, baseReactant and baseProduct in cellDesigner
#R2 <- readLines(file("txt6_repeated_part2")) # define other repeated Reactants in cellDesigner
R3 <- readLines(file("txt6_repeated_part3")) # connect two parts
#R4 <- readLines(file("txt6_repeated_part4")) # define other repeated products in cellDesigner
#R5 <- readLines(file("txt6_repeated_part5")) # mainly define connectScheme
R5 <- readLines(file("txt6_repeated_part5.xml")) # mainly define connectScheme____version2
#R6 <- readLines(file("txt6_repeated_part6")) # define the repeated reactants 
R7 <- readLines(file("txt6_repeated_part7")) # connect two parts
#R8 <- readLines(file("txt6_repeated_part8")) # define the repeated products
#R9 <- readLines(file("txt6_repeated_part9")) # the end of reaction annotation
#R9 <- readLines(file("txt6_repeated_part9.xml")) # the end of reaction annotation___version2

# the output of the followed function could be replaced of R1
getReactionMain <- function(rxnID, BR_specie, BR_id, BP_specie, BP_id ){
  txt0 <- vector()
  txt0[1]  <- paste("<reaction metaid=\"", rxnID, "\" id=\"", rxnID, "\" reversible=\"false\">", sep = "")             
  txt0[2]  <- "<annotation>"                                                           
  txt0[3]  <- "<celldesigner:extension>"                                               
  txt0[4]  <- "<celldesigner:reactionType>STATE_TRANSITION</celldesigner:reactionType>"
  txt0[5]  <- "<celldesigner:baseReactants>"                                           
  txt0[6]  <- paste("<celldesigner:baseReactant species=\"", BR_specie, "\" alias=\"", BR_id, "\">", sep = "")            
  txt0[7]  <- "<celldesigner:linkAnchor position=\"E\"/>"                              
  txt0[8]  <- "</celldesigner:baseReactant>"                                           
  txt0[9]  <- "</celldesigner:baseReactants>"                                          
  txt0[10]  <- "<celldesigner:baseProducts>"                                            
  txt0[11]  <- paste("<celldesigner:baseProduct species=\"",BP_specie, "\" alias=\"", BP_id, "\">", sep = "")               
  txt0[12]  <- "<celldesigner:linkAnchor position=\"W\"/>"                              
  txt0[13]  <- "</celldesigner:baseProduct>"                                            
  txt0[14]  <- "</celldesigner:baseProducts>"                                           
  txt0[15]  <- "<celldesigner:listOfReactantLinks>"
  return(txt0)
  
}

# the output of the followed function could be replaced of R2
getOtherReactant <- function(R_specie, R_id ){
  txt0 <- vector()
  txt0[1] <- paste("<celldesigner:reactantLink reactant=\"",R_specie,"\" alias=\"",R_id,"\" targetLineIndex=\"-1,0\">", sep="")
  txt0[2] <- "<celldesigner:linkAnchor position=\"E\"/>"                                          
  txt0[3] <- "<celldesigner:connectScheme connectPolicy=\"direct\">"                              
  txt0[4] <- "<celldesigner:listOfLineDirection>"                                                 
  txt0[5] <- "<celldesigner:lineDirection index=\"0\" value=\"unknown\"/>"                        
  txt0[6] <- "</celldesigner:listOfLineDirection>"                                                
  txt0[7] <- "</celldesigner:connectScheme>"                                                      
  txt0[8] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\" type=\"Straight\"/>"            
  txt0[9] <- "</celldesigner:reactantLink>" 
  return(txt0)
}

getOtherProduct <- function(P_specie, P_id ){
  txt0 <- vector()
  txt0[1] <- paste("<celldesigner:productLink product=\"",P_specie,"\" alias=\"",P_id,"\" targetLineIndex=\"-1,1\">",sep="")
  txt0[2] <- "<celldesigner:linkAnchor position=\"W\"/>"                                        
  txt0[3] <- "<celldesigner:connectScheme connectPolicy=\"direct\">"                            
  txt0[4] <- "<celldesigner:listOfLineDirection>"                                               
  txt0[5] <- "<celldesigner:lineDirection index=\"0\" value=\"unknown\"/>"                      
  txt0[6] <- "</celldesigner:listOfLineDirection>"                                              
  txt0[7] <- "</celldesigner:connectScheme>"                                                    
  txt0[8] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\" type=\"Straight\"/>"          
  txt0[9] <- "</celldesigner:productLink>"
  
  return(txt0)
}

DefineAllMet <- function(specie, id, Metaid){
  txt0 <- vector()
  txt0[1] <- paste("<speciesReference metaid=\"",Metaid,"\" species=\"",specie,"\">", sep="")
  txt0[2] <- "<annotation>"                                           
  txt0[3] <- "<celldesigner:extension>"                               
  txt0[4] <- paste("<celldesigner:alias>",id,"</celldesigner:alias>", sep="")          
  txt0[5] <- "</celldesigner:extension>"                              
  txt0[6] <- "</annotation>"                                          
  txt0[7] <- "</speciesReference>"  
  return(txt0)
}

DeProteinConnect <- function(p_specie, p_id, g_specie, g_id){
  txt0 <- vector()
  txt0[1] <- "</celldesigner:listOfProductLinks>"                                                                        
  txt0[2] <- "<celldesigner:connectScheme connectPolicy=\"direct\" rectangleIndex=\"0\">"                                
  txt0[3] <- "<celldesigner:listOfLineDirection>"                                                                        
  txt0[4] <- "<celldesigner:lineDirection index=\"0\" value=\"unknown\"/>"                                               
  txt0[5] <- "</celldesigner:listOfLineDirection>"                                                                       
  txt0[6] <- "</celldesigner:connectScheme>"                                                                             
  txt0[7] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\"/>"                                                     
  txt0[8] <- "<celldesigner:listOfModification>"                                                                         
  txt0[9] <- paste("<celldesigner:modification type=\"CATALYSIS\" modifiers=\"",p_specie,"\" aliases=\"",p_id,"\" targetLineIndex=\"-1,4\">", sep = "")
  txt0[10] <- "<celldesigner:connectScheme connectPolicy=\"direct\">"                                                     
  txt0[11] <- "<celldesigner:listOfLineDirection>"                                                                        
  txt0[12] <- "<celldesigner:lineDirection index=\"0\" value=\"unknown\"/>"                                               
  txt0[13] <- "</celldesigner:listOfLineDirection>"                                                                       
  txt0[14] <- "</celldesigner:connectScheme>"                                                                             
  txt0[15] <- paste("<celldesigner:linkTarget species=\"",p_specie,"\" alias=\"",p_id,"\">", sep = "")                                                  
  txt0[16] <- "<celldesigner:linkAnchor position=\"S\"/>"                                                                 
  txt0[17] <- "</celldesigner:linkTarget>"                                                                                
  txt0[18] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\"/>"                                                     
  txt0[19] <- "</celldesigner:modification>"                                                                              
  txt0[20] <- paste("<celldesigner:modification type=\"CATALYSIS\" modifiers=\"",g_specie,"\" aliases=\"",g_id,"\" targetLineIndex=\"-1,5\">", sep = "")
  txt0[21] <- "<celldesigner:connectScheme connectPolicy=\"direct\">"                                                     
  txt0[22] <- "<celldesigner:listOfLineDirection>"                                                                        
  txt0[23] <- "<celldesigner:lineDirection index=\"0\" value=\"unknown\"/>"                                               
  txt0[24] <- "</celldesigner:listOfLineDirection>"                                                                       
  txt0[25] <- "</celldesigner:connectScheme>"                                                                             
  txt0[26] <- paste("<celldesigner:linkTarget species=\"",g_specie,"\" alias=\"",g_id,"\">", sep = "")                                
  txt0[27] <- "<celldesigner:linkAnchor position=\"S\"/>"                                                                 
  txt0[28] <- "</celldesigner:linkTarget>"                                                                                
  txt0[29] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\"/>"                                                     
  txt0[30] <- "</celldesigner:modification>"                                                                              
  txt0[31] <- "</celldesigner:listOfModification>"                                                                        
  txt0[32] <- "</celldesigner:extension>"                                                                                 
  txt0[33] <- "</annotation>"                                                                                             
  txt0[34] <- "<listOfReactants>" 
  return(txt0)
}


DefineProtein <- function(p_specie, p_id, p_metaid, g_specie, g_id, g_metaid){
  txt0 <- vector()
  txt0[1] <- "</listOfProducts>"                                              
  txt0[2] <- "<listOfModifiers>"                                              
  txt0[3] <- paste("<modifierSpeciesReference metaid=\"",p_metaid,"\" species=\"",p_specie,"\">", sep = "")
  txt0[4] <- "<annotation>"                                                   
  txt0[5] <- "<celldesigner:extension>"                                       
  txt0[6] <- paste("<celldesigner:alias>",p_id,"</celldesigner:alias>", sep = "")                 
  txt0[7] <- "</celldesigner:extension>"                                      
  txt0[8] <- "</annotation>"                                                  
  txt0[9] <- "</modifierSpeciesReference>"                                    
  txt0[10] <- paste("<modifierSpeciesReference metaid=\"",g_metaid,"\" species=\"",g_specie,"\">", sep = "")
  txt0[11] <- "<annotation>"                                                   
  txt0[12] <- "<celldesigner:extension>"                                       
  txt0[13] <- paste("<celldesigner:alias>",g_id,"</celldesigner:alias>", sep = "")                 
  txt0[14] <- "</celldesigner:extension>"                                      
  txt0[15] <- "</annotation>"                                                  
  txt0[16] <- "</modifierSpeciesReference>"                                    
  txt0[17] <- "</listOfModifiers>"                                             
  txt0[18] <- "</reaction>"
  return(txt0)
}



getRxnInformation <- function (rxn0, rxn_p, R30=R3,R70=R7){ #input a datafram contains the detailed reaction information
  #define the reaction metaid, baseReactant and baseProduct in cellDesigner
  #test
  #rxn0 = rxn_annotation[[1]]
  #rxn_p = gpr[[1]]
  #rxn0$note[3]
  
  rxnID <- rxn0$rxnID[1]
  BR_specie <- rxn0$specie[ which(rxn0$type=="reactant" & rxn0$note=="base")]
  BR_id <- rxn0$id[ which(rxn0$type=="reactant" & rxn0$note=="base")]
  BP_specie <- rxn0$specie[ which(rxn0$type=="product" & rxn0$note=="base")]
  BP_id <- rxn0$id[ which(rxn0$type=="product" & rxn0$note=="base")]
  txt6_R1 <- getReactionMain(rxnID,BR_specie,BR_id,BP_specie,BP_id)
  
  #define other repeated Reactants in cellDesigner
  if(length(which(rxn0$type=="reactant" & rxn0$note != "base"))){
    otherReactantID <- which(rxn0$type=="reactant" & rxn0$note != "base") ##should be careful about this result
    txt6_R2 <- vector()
    R_specie <- rxn0$specie[otherReactantID]
    R_id <- rxn0$id[otherReactantID]
    txt6_R2 <- vector()
    for (i in seq(length((otherReactantID)))){
      txt6_R2 <- c(txt6_R2,getOtherReactant(R_specie[i], R_id[i]))
    }
  } else{
    txt6_R2 <- NULL
  }
  
  txt6_R3 <- R30
  #define other repeated Products in cellDesigner  
  if(length(which(rxn0$type=="product" & rxn0$note != "base"))){
    otherProductID <- which(rxn0$type=="product" & rxn0$note != "base")
    txt6_R4 <- vector()
    P_specie <- rxn0$specie[otherProductID]
    P_id <- rxn0$id[otherProductID]
    txt6_R4 <- vector()
    for (i in seq(length(otherProductID))){
      txt6_R4 <- c(txt6_R4,getOtherProduct(P_specie[i], P_id[i]))
    }
  } else{
    txt6_R4 <- NULL
  }
  

  txt6_R5 <- DeProteinConnect(rxn_p$protein_specie,rxn_p$protein_id, rxn_p$gene_specie,rxn_p$gene_id)
  
  #define all reactants information
  if(length(which(rxn0$type=="reactant"))){
    reactantID <- which(rxn0$type=="reactant")
    specie_r <- rxn0$specie[reactantID]
    id_r <- rxn0$id[reactantID]
    MetaID_r <- rxn0$MetaID[reactantID]
    txt6_R6 <- vector()
    for (i in seq(length(reactantID))){
      txt6_R6 <- c(txt6_R6,DefineAllMet(specie_r[i], id_r[i], MetaID_r[i]))
    }
  } else{
    txt6_R6 <- ""
  }
  
  txt6_R7 <-R70
  
  #define all products information
  if(length(which(rxn0$type=="product"))){
    productID <- which(rxn0$type=="product")
    specie_p <- rxn0$specie[productID]
    id_p <- rxn0$id[productID]
    MetaID_p <- rxn0$MetaID[productID]
    txt6_R8 <- vector()
    for (i in seq(length(productID))){
      txt6_R8 <- c(txt6_R8,DefineAllMet(specie_p[i], id_p[i], MetaID_p[i]))
    }
  } else{
    txt6_R8 <- ""
  }
  
  txt6_R9 <- DefineProtein(rxn_p$protein_specie,rxn_p$protein_id, rxn_p$MetaID_p, rxn_p$gene_specie,rxn_p$gene_id, rxn_p$MetaID_g)
  
  
  txt6_new <- c(txt6_R1,txt6_R2,txt6_R3,txt6_R4,txt6_R5,txt6_R6,txt6_R7,txt6_R8,txt6_R9)

  return(txt6_new)
  
  }



## test 1  batch process
##import all reaction
rxn <- read_excel("rxnInformation_add_reaction_automatic.xlsx")
rxn_p0 <- read_excel("rxnInformation_add_reaction_automatic.xlsx", sheet = "gpr")

##change reaction into list format
## data process before using function
rxnID <- unique(rxn$rxnID)
rxn_annotation <- list()
rxn_id <- list()
for (i in 1:length(rxnID)){
  rxn_id[[i]] <- which(rxn$rxnID %in% rxnID[i] == TRUE)
  rxn_annotation[[i]] <- rxn[rxn_id[[i]],]
  }

gpr <- list()
for (i in 1:length(rxn_p0$rxnID)){
gpr[[i]] <- rxn_p0[i,]
}

## using the main function
txt6_new <- vector()
for (j in 1:length(rxnID)){
txt6_new <- c(txt6_new,getRxnInformation (rxn_annotation[[j]],rxn_p = gpr[[j]]))
}


txt_new <- c(txt1,txt2,txt30,txt4,txt5,txt6_new,txt7) # connect each template
writeLines(txt_new, file("template_test_protein-add rxn connection_using function.xml")) # save





###application 1
rxn <- rxn_core_carbon_cellD0
rxn_p0 <- gpr

##change reaction into list format

## data process before using function
rxnID <- unique(rxn$rxnID)
rxn_annotation <- list()
rxn_id <- list()
for (i in 1:length(rxnID)){
  rxn_id[[i]] <- which(rxn$rxnID %in% rxnID[i] == TRUE)
  rxn_annotation[[i]] <- rxn[rxn_id[[i]],]
}

gpr <- list()
for (i in 1:length(rxn_p0$rxnID)){
  gpr[[i]] <- rxn_p0[i,]
}

## using the main function
txt6_new <- vector()
for (j in 1:length(rxnID)){
  txt6_new <- c(txt6_new,getRxnInformation (rxn_annotation[[j]],rxn_p = gpr[[j]]))
}

txt_new <- c(txt1,txt2,txt30,txt4,txt5,txt6_new,txt7) # connect each template
writeLines(txt_new, file("model_test_protein-add rxn connection_using function.xml")) # save










