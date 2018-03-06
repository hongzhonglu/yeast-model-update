##txt6 is the main reaction format information from xml of cellDesigner
##txt6 <- readLines(file("template6")) # repeated part: information of each reactions
##split txt6_repeated_part into several part
#R1 <- readLines(file("txt6_repeated_part1")) # define the reaction metaid, baseReactant and baseProduct in cellDesigner
#R2 <- readLines(file("txt6_repeated_part2")) # define other repeated Reactants in cellDesigner
R3 <- readLines(file("txt6_repeated_part3")) # connect two parts
#R4 <- readLines(file("txt6_repeated_part4")) # define other repeated products in cellDesigner
R5 <- readLines(file("txt6_repeated_part5")) # mainly define connectScheme
#R6 <- readLines(file("txt6_repeated_part6")) # define the repeated reactants 
R7 <- readLines(file("txt6_repeated_part7")) # connect two parts
#R8 <- readLines(file("txt6_repeated_part8")) # define the repeated products
R9 <- readLines(file("txt6_repeated_part9")) # the end of reaction annotation

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


getRxnInformation <- function (rxn0, R30=R3,R50=R5,R70=R7,R90=R9){ #input a datafram contains the detailed reaction information
  #define the reaction metaid, baseReactant and baseProduct in cellDesigner
  rxnID <- rxn0$rxnID[1]
  BR_specie <- rxn0$specie[ which(rxn0$type=="reactant" & rxn0$note=="base")]
  BR_id <- rxn0$id[ which(rxn0$type=="reactant" & rxn0$note=="base")]
  BP_specie <- rxn0$specie[ which(rxn0$type=="product" & rxn0$note=="base")]
  BP_id <- rxn0$id[ which(rxn0$type=="product" & rxn0$note=="base")]
  txt6_R1 <- getReactionMain(rxnID,BR_specie,BR_id,BP_specie,BP_id)
  
  #define other repeated Reactants in cellDesigner
  if(length(which(rxn0$type=="reactant" & is.na(rxn0$note)))){
    otherReactantID <- which(rxn0$type=="reactant" & is.na(rxn0$note))
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
  if(length(which(rxn0$type=="product" & is.na(rxn0$note)))){
    otherProductID <- which(rxn0$type=="product" & is.na(rxn0$note))
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
  
  txt6_R5 <-R50
  
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
  
  txt6_R9 <- R90
  txt6_new <- c(txt6_R1,txt6_R2,txt6_R3,txt6_R4,txt6_R5,txt6_R6,txt6_R7,txt6_R8,txt6_R9)
}
## use function to get the txt6--the repeated part of the reactions in cell designer





## test 1  batch process
##import all reaction
rxn <- read_excel("rxnInformation_add_reaction_automatic.xlsx")
##change reaction into list format
rxnID <- unique(rxn$rxnID)
rxn_annotation <- list()
rxn_id <- list()
for (i in 1:length(rxnID)){
  rxn_id[[i]] <- which(rxn$rxnID %in% rxnID[i] == TRUE)
  rxn_annotation[[i]] <- rxn[rxn_id[[i]],]
  }

txt6_new <- vector()
for (j in 1:length(rxnID)){
txt6_new <- c(txt6_new,getRxnInformation (rxn_annotation[[j]]))
}

txt_new <- c(txt1,txt2,txt3,txt4,txt5,txt6_new,txt7) # connect each template
writeLines(txt_new, file("template_test_add rxn connection_using function.xml")) # save



###application 1 import reactions from transport subsystems
##import all reaction
rxn <- read_excel("rxn_from_transport_yeastGEM.xlsx")
rxn <- rxn[1:336,]
##change reaction into list format
rxnID <- unique(rxn$rxnID)
rxn_annotation <- list()
rxn_id <- list()
for (i in 1:length(rxnID)){
  rxn_id[[i]] <- which(rxn$rxnID %in% rxnID[i] == TRUE)
  rxn_annotation[[i]] <- rxn[rxn_id[[i]],]
}

txt6_new <- vector()
for (j in 1:length(rxnID)){
  txt6_new <- c(txt6_new,getRxnInformation (rxn_annotation[[j]]))
}

txt_new <- c(txt1,txt2,txt3,txt4,txt5,txt6_new,txt7) # connect each template
writeLines(txt_new, file("model_add rxn connection_using function.xml")) # save



###application 2 import reactions from core carbon metabolism: EMP,PPP and TCA cycles
##import all reaction
rxn <- read_excel("rxn_from_core_carbon_yeastGEM.xlsx")


#remove currency
which(str_detect(rxn$name,"H\\+")!=TRUE)
rxn <- rxn[which(str_detect(rxn$name,"H\\+")!=TRUE),]
rxn <- rxn[which(str_detect(rxn$name,"H2O")!=TRUE),]


##change reaction into list format
rxnID <- unique(rxn$rxnID)
rxn_annotation <- list()
rxn_id <- list()
for (i in 1:length(rxnID)){
  rxn_id[[i]] <- which(rxn$rxnID %in% rxnID[i] == TRUE)
  rxn_annotation[[i]] <- rxn[rxn_id[[i]],]
}

txt6_new <- vector()
for (j in 1:length(rxnID)){
  txt6_new <- c(txt6_new,getRxnInformation (rxn_annotation[[j]]))
}

txt_new <- c(txt1,txt2,txt3,txt4,txt5,txt6_new,txt7) # connect each template
writeLines(txt_new, file("model_add rxn connection_using function2.xml")) # save
