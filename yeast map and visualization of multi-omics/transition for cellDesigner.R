#-----------------------------------------------------------------------
#prepare the component for the cellDesigner
#------------------------------------------------------------------------

# define the size of graph
getNewSize <- function(x=600,y=400){
  txt <- vector()
  txt[1] <-  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"                                                                                                                                                                                                                
  txt[2] <-  "<sbml xmlns=\"http://www.sbml.org/sbml/level2/version4\" xmlns:celldesigner=\"http://www.sbml.org/2001/ns/celldesigner\" level=\"2\" version=\"4\">"                                                                                                       
  txt[3] <-  "<model metaid=\"Acetate_space_utilization\" id=\"Acetate_space_utilization\">"                                                                                                                                                                             
  txt[4] <-  "<notes>"                                                                                                                                                                                                                                                   
  txt[5] <-  "<html xmlns=\"http://www.w3.org/1999/xhtml\">"                                                                                                                                                                                                             
  txt[6] <-  "<head>"                                                                                                                                                                                                                                                    
  txt[7] <-  "<title/>"                                                                                                                                                                                                                                                  
  txt[8] <-  "</head>"                                                                                                                                                                                                                                                   
  txt[9] <-  "<body>The prokaryotic utilization of acetate and its interconversion to acetyl-CoA.  Acetate and acetyl-CoA can function as a carbon source for energy production. Website=www.biocyc.org/Ecoli/Pathways/free text=Voet and Voet  (textbook)  Biochemistry"
  txt[10] <-  "</body>"                                                                                                                                                                                                                                                   
  txt[11] <-  "</html>"                                                                                                                                                                                                                                                   
  txt[12] <-  "</notes>"                                                                                                                                                                                                                                                  
  txt[13] <-  "<annotation>"                                                                                                                                                                                                                                              
  txt[14] <-  "<celldesigner:extension>"                                                                                                                                                                                                                                  
  txt[15] <-  "<celldesigner:modelVersion>4.0</celldesigner:modelVersion>"                                                                                                                                                                                                
  txt[16] <-  paste("<celldesigner:modelDisplay sizeX=\"",x,"\" sizeY=\"",y,"\"/>",sep="")                                                                                                                                                                                                  
  txt[17] <-  "<celldesigner:listOfCompartmentAliases/>"                                                                                                                                                                                                                  
  txt[18] <-  "<celldesigner:listOfComplexSpeciesAliases/>"                                                                                                                                                                                                               
  txt[19] <-  "<celldesigner:listOfSpeciesAliases>" 
  return(txt)
}

# define position information for each metabolite
getPosition <- function(id01,species01,x01,y01){
  #input
  #id01 id of a component
  #species01 species id of a component
  #x01 x coordinate of a component
  #y01 y coordinate of a component
  #output
  #a text file with coordinate of the component in a map
  txt <- vector()
  txt[1] <- paste("<celldesigner:speciesAlias id=\"",id01,"\" species=\"",species01,"\">", sep = "")           
  txt[2] <- "<celldesigner:activity>inactive</celldesigner:activity>"            
  txt[3] <- paste("<celldesigner:bounds x=\"",x01,"\" y=\"",y01,"\" w=\"70.0\" h=\"25.0\"/>", sep = "")
  txt[4] <- "<celldesigner:font size=\"12\"/>"                                   
  txt[5] <- "<celldesigner:view state=\"usual\"/>"                               
  txt[6] <- "<celldesigner:usualView>"                                           
  txt[7] <- "<celldesigner:innerPosition x=\"0.0\" y=\"0.0\"/>"                  
  txt[8] <- "<celldesigner:boxSize width=\"70.0\" height=\"25.0\"/>"             
  txt[9] <- "<celldesigner:singleLine width=\"1.0\"/>"                           
  txt[10] <- "<celldesigner:paint color=\"ffccff66\" scheme=\"Color\"/>"          
  txt[11] <- "</celldesigner:usualView>"                                          
  txt[12] <- "<celldesigner:briefView>"                                           
  txt[13] <- "<celldesigner:innerPosition x=\"0.0\" y=\"0.0\"/>"                  
  txt[14] <- "<celldesigner:boxSize width=\"80.0\" height=\"60.0\"/>"             
  txt[15] <- "<celldesigner:singleLine width=\"0.0\"/>"                           
  txt[16] <- "<celldesigner:paint color=\"3fff0000\" scheme=\"Color\"/>"          
  txt[17] <- "</celldesigner:briefView>"                                          
  txt[18] <- "<celldesigner:info state=\"empty\" angle=\"-1.5707963267948966\"/>" 
  txt[19] <- "</celldesigner:speciesAlias>" 
  return(txt)
}


# define annotation information for each metabolite
getAnnotation <- function(species01,name01){
  #input
  #species01 species id of a component
  #name01 name of a component
  #output
  #a text file with name of the component in a map
  
  ##class_type: SIMPLE_MOLECULE,PROTEIN, GENE
  txt <- vector()
  txt[1] <- paste("<species metaid=\"",species01,"\" id=\"",species01,"\" name=\"",name01,"\" compartment=\"default\" initialAmount=\"0\">",sep = "")             
  txt[2] <- "<annotation>"                                                                                            
  txt[3] <- "<celldesigner:extension>"                                                                                
  txt[4] <- "<celldesigner:positionToCompartment>inside</celldesigner:positionToCompartment>"                         
  txt[5] <- "<celldesigner:speciesIdentity>"                                                                          
  txt[6] <- "<celldesigner:class>SIMPLE_MOLECULE</celldesigner:class>"                                            
  txt[7] <- paste("<celldesigner:name>",name01,"</celldesigner:name>", sep = "")                                                             
  txt[8] <- "</celldesigner:speciesIdentity>"                                                                         
  txt[9] <- "</celldesigner:extension>"                                                                               
  txt[10] <- "</annotation>"                                                                                           
  txt[11] <- "</species>"  
  return(txt)
}

# define annotation information for each protein
getProAnnotation <- function(species01,name01, proteinID01){
  #input
  #species01 species id of a component
  #name01 name of a component
  #proteinID01 ID of the protein as the component
  #output
  #a text file with protein annotation for the related component
  
  
  ##class_type: SIMPLE_MOLECULE,PROTEIN, GENE
  txt <- vector()
  txt[1] <- paste("<species metaid=\"",species01,"\" id=\"",species01,"\" name=\"",name01,"\" compartment=\"default\" initialAmount=\"0\">", sep="")
  txt[2] <- "<annotation>"                                                                             
  txt[3] <- "<celldesigner:extension>"                                                                 
  txt[4] <- "<celldesigner:positionToCompartment>inside</celldesigner:positionToCompartment>"          
  txt[5] <- "<celldesigner:speciesIdentity>"                                                           
  txt[6] <- "<celldesigner:class>PROTEIN</celldesigner:class>"                                         
  txt[7] <- paste("<celldesigner:proteinReference>",proteinID01,"</celldesigner:proteinReference>", sep="")                       
  txt[8] <- "</celldesigner:speciesIdentity>"                                                          
  txt[9] <- "</celldesigner:extension>"                                                                
  txt[10] <- "</annotation>"                                                                            
  txt[11] <- "</species>"
  return(txt)
}

# define annotation information for each gene
getGeneAnnotation <- function(species01,name01, geneID01){
  #input
  #species01 species id of a component
  #name01 name of a component
  #geneID01 ID of the gene as the component
  #output
  #a text file with gene annotation for the related component
  
  txt <- vector()
  txt[1] <- paste("<species metaid=\"",species01,"\" id=\"",species01,"\" name=\"",name01,"\" compartment=\"default\" initialAmount=\"0\">", sep = "")
  txt[2] <- "<annotation>"                                                                             
  txt[3] <- "<celldesigner:extension>"                                                                 
  txt[4] <- "<celldesigner:positionToCompartment>inside</celldesigner:positionToCompartment>"          
  txt[5] <- "<celldesigner:speciesIdentity>"                                                           
  txt[6] <- "<celldesigner:class>GENE</celldesigner:class>"                                            
  txt[7] <- paste("<celldesigner:geneReference>",geneID01,"</celldesigner:geneReference>", sep="")                             
  txt[8] <- "</celldesigner:speciesIdentity>"                                                          
  txt[9] <- "</celldesigner:extension>"                                                                
  txt[10] <- "</annotation>"                                                                            
  txt[11] <- "</species>" 
  return(txt)
}


#-----------------------------------------------------------------------
#prepare the rxn connection for the cellDesigner
#------------------------------------------------------------------------

# the output of the followed function could be replaced of R1
getReactionMain <- function(rxnID, BR_specie, BR_id, BP_specie, BP_id ){
  txt0 <- vector()
  txt0[1]  <- paste("<reaction metaid=\"", rxnID, "\" id=\"", rxnID, "\" reversible=\"false\">", sep = "")             
  txt0[2]  <- "<annotation>"                                                           
  txt0[3]  <- "<celldesigner:extension>"                                               
  txt0[4]  <- "<celldesigner:reactionType>STATE_TRANSITION</celldesigner:reactionType>"
  txt0[5]  <- "<celldesigner:baseReactants>"                                           
  txt0[6]  <- paste("<celldesigner:baseReactant species=\"", BR_specie, "\" alias=\"", BR_id, "\">", sep = "")            
  #txt0[7]  <- "<celldesigner:linkAnchor position=\"E\"/>"     ##define how the metabolite is connected??                         
  txt0[7]  <- "</celldesigner:baseReactant>"                                           
  txt0[8]  <- "</celldesigner:baseReactants>"                                          
  txt0[9]  <- "<celldesigner:baseProducts>"                                            
  txt0[10]  <- paste("<celldesigner:baseProduct species=\"",BP_specie, "\" alias=\"", BP_id, "\">", sep = "")               
  #txt0[11]  <- "<celldesigner:linkAnchor position=\"W\"/>"     ##define how the metabolite is connected??                             
  txt0[11]  <- "</celldesigner:baseProduct>"                                            
  txt0[12]  <- "</celldesigner:baseProducts>"                                           
  txt0[13]  <- "<celldesigner:listOfReactantLinks>"
  return(txt0)
  
}

# the output of the followed function could be replaced of R2
getOtherReactant <- function(R_specie, R_id ){
  txt0 <- vector()
  txt0[1] <- paste("<celldesigner:reactantLink reactant=\"",R_specie,"\" alias=\"",R_id,"\" targetLineIndex=\"-1,0\">", sep="")
  #txt0[2] <- "<celldesigner:linkAnchor position=\"E\"/>"    ##define how the metabolite is connected??                                                 
  txt0[2] <- "<celldesigner:connectScheme connectPolicy=\"direct\">"                              
  txt0[3] <- "<celldesigner:listOfLineDirection>"                                                 
  txt0[4] <- "<celldesigner:lineDirection index=\"0\" value=\"unknown\"/>"                        
  txt0[5] <- "</celldesigner:listOfLineDirection>"                                                
  txt0[6] <- "</celldesigner:connectScheme>"                                                      
  txt0[7] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\" type=\"Straight\"/>"            
  txt0[8] <- "</celldesigner:reactantLink>" 
  return(txt0)
}

getOtherProduct <- function(P_specie, P_id ){
  txt0 <- vector()
  txt0[1] <- paste("<celldesigner:productLink product=\"",P_specie,"\" alias=\"",P_id,"\" targetLineIndex=\"-1,1\">",sep="")
  #txt0[2] <- "<celldesigner:linkAnchor position=\"W\"/>"     ##define how the metabolite is connected??                                      
  txt0[2] <- "<celldesigner:connectScheme connectPolicy=\"direct\">"                            
  txt0[3] <- "<celldesigner:listOfLineDirection>"                                               
  txt0[4] <- "<celldesigner:lineDirection index=\"0\" value=\"unknown\"/>"                      
  txt0[5] <- "</celldesigner:listOfLineDirection>"                                              
  txt0[6] <- "</celldesigner:connectScheme>"                                                    
  txt0[7] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\" type=\"Straight\"/>"          
  txt0[8] <- "</celldesigner:productLink>"
  
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
  #txt0[16] <- "<celldesigner:linkAnchor position=\"S\"/>"   ##define how the metabolite is connected??                                                                 
  txt0[16] <- "</celldesigner:linkTarget>"                                                                                
  txt0[17] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\"/>"                                                     
  txt0[18] <- "</celldesigner:modification>"                                                                              
  txt0[19] <- paste("<celldesigner:modification type=\"CATALYSIS\" modifiers=\"",g_specie,"\" aliases=\"",g_id,"\" targetLineIndex=\"-1,5\">", sep = "")
  txt0[20] <- "<celldesigner:connectScheme connectPolicy=\"direct\">"                                                     
  txt0[21] <- "<celldesigner:listOfLineDirection>"                                                                        
  txt0[22] <- "<celldesigner:lineDirection index=\"0\" value=\"unknown\"/>"                                               
  txt0[23] <- "</celldesigner:listOfLineDirection>"                                                                       
  txt0[24] <- "</celldesigner:connectScheme>"                                                                             
  txt0[25] <- paste("<celldesigner:linkTarget species=\"",g_specie,"\" alias=\"",g_id,"\">", sep = "")                                
  #txt0[27] <- "<celldesigner:linkAnchor position=\"S\"/>"   ##define how the metabolite is connected??                                                               
  txt0[26] <- "</celldesigner:linkTarget>"                                                                                
  txt0[27] <- "<celldesigner:line width=\"1.0\" color=\"ff000000\"/>"                                                     
  txt0[28] <- "</celldesigner:modification>"                                                                              
  txt0[29] <- "</celldesigner:listOfModification>"                                                                        
  txt0[30] <- "</celldesigner:extension>"                                                                                 
  txt0[31] <- "</annotation>"                                                                                             
  txt0[32] <- "<listOfReactants>" 
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


produceInputForCellDesigner <- function(met_annotation_inf, 
                                        gpr_inf, 
                                        rxn_core_carbon_inf,
                                        x_size=2500, 
                                        y_size=3600) {
  # this script is used to input the coordinate of metaoblite, protein and gene into celldesigner
  # Input the comon template
  txt1 <- readLines(file("data/template1")) # model information and size
  txt3 <- readLines(file("data/template3")) # other information of model besides metabolits and reactions
  txt5 <- readLines(file("data/template5")) # sign of xml format
  txt7 <- readLines(file("data/template7")) # end sign of xml format
  
  txt1 <- getNewSize(x = 2500, y = 3600)
  txt2 <- vector()
  txt3 <- vector()
  txt4 <- vector()
  
  #--------------------------------------
  # txt2
  #--------------------------------------
  id <- met_annotation_inf$id
  species <- met_annotation_inf$species
  name <- met_annotation_inf$name
  x <- as.character(met_annotation_inf$x)
  y <- as.character(met_annotation_inf$y)
  
  for (i in seq(length(id))) {
    txt2 <- c(getPosition(id01 = id[i], species01 = species[i], x01 = x[i], y01 = y[i]), txt2)
  }
  
  #--------------------------------------
  # txt3
  #--------------------------------------
  ## add protein define
  ## add gene define
  ## input protein and gene information
  protein0 <- filter(met_annotation_inf, type == "PROTEIN")
  protein0$type <- "GENERIC"
  # protein0$protein_id <- paste("pr",1:length(protein0$metID), sep="")
  protein0$protein_id <- str_replace_all(protein0$rxnID, "r_", "pr")
  
  
  gene0 <- filter(met_annotation_inf, type == "GENE")
  # gene0$gene_id <- paste("gn",1:length(gene0$metID), sep="")
  gene0$gene_id <- str_replace_all(gene0$rxnID, "r_", "gn")
  
  getProteinGene <- function(id, type, name) {
    txt0 <- vector()
    if (type == "GENERIC") {
      txt0[1] <- paste("<celldesigner:protein id=\"", id, "\" name=\"", name, "\" type=\"GENERIC\"/>", sep = "")
    }
    if (type == "GENE") {
      txt0[1] <- paste("<celldesigner:gene id=\"", id, "\" name=\"", name, "\" type=\"GENE\"/>", sep = "")
    }
    
    return(txt0)
  }
  
  
  txt3_1 <- readLines(file("data/txt3_1.xml"))
  
  # txt3_2 protein module
  protein_id <- protein0$protein_id
  type_p <- protein0$type
  name_p <- protein0$name
  
  
  txt3_2 <- vector()
  for (i in seq(length(protein_id))) {
    txt3_2 <- c(txt3_2, getProteinGene(id = protein_id[i], type = type_p[i], name = name_p[i]))
  }
  
  
  txt3_3 <- readLines(file("data/txt3_3.xml"))
  
  # txt3_4 gene module
  gene_id <- gene0$gene_id
  type_g <- gene0$type
  name_g <- gene0$name
  txt3_4 <- vector()
  for (i in seq(length(gene_id))) {
    txt3_4 <- c(txt3_4, getProteinGene(id = gene_id[i], type = type_g[i], name = name_g[i]))
  }
  
  txt3_5 <- readLines(file("data/txt3_5.xml"))
  txt30 <- c(txt3_1, txt3_2, txt3_3, txt3_4, txt3_5)
  
  
  #--------------------------------------------------------------------------
  ## txt4: give the annotation for each unique metabolite
  #--------------------------------------------------------------------------
  # txt4_1 metabolite module
  met_without_protein0 <- filter(met_annotation_inf, type == "SIMPLE_MOLECULE")
  species_unique <- unique(met_without_protein0$species)
  name_unique <- unique(met_without_protein0$name)
  txt4_1 <- vector()
  for (i in seq(length(species_unique))) {
    txt4_1 <- c(txt4_1, getAnnotation(species01 = species_unique[i], name01 = name_unique[i]))
  }
  
  # txt4_2 protein module
  protein_specie <- protein0$species
  protein_id <- protein0$protein_id
  type_p <- protein0$type
  name_p <- protein0$name
  
  txt4_2 <- vector()
  for (i in seq(length(protein_id))) {
    txt4_2 <- c(txt4_2, getProAnnotation(species01 = protein_specie[i], name01 = name_p[i], proteinID01 = protein_id[i]))
  }
  
  # txt4_3 gene module
  gene_id <- gene0$gene_id
  gene_specie <- gene0$species
  type_g <- gene0$type
  name_g <- gene0$name
  txt4_3 <- vector()
  for (i in seq(length(gene_id))) {
    txt4_3 <- c(txt4_3, getGeneAnnotation(species01 = gene_specie[i], name01 = name_g[i], geneID01 = gene_id[i]))
  }
  txt4 <- c(txt4_1, txt4_2, txt4_3)
  
  
  # merge the above file
  txt_final0 <- c(txt1, txt2, txt30, txt4, txt5, txt7)
  
  #----------------------------------------------------------
  # txt6 Input the connection of metabolite into celldesigner
  #---------------------------------------------------------
  # input the common template
  R3 <- readLines(file("data/txt6_repeated_part3")) # connect two parts
  R5 <- readLines(file("data/txt6_repeated_part5.xml")) # mainly define connectScheme____version2
  R7 <- readLines(file("data/txt6_repeated_part7")) # connect two parts
  
  rxn <- rxn_core_carbon_inf
  rxn_p0 <- gpr_inf
  ## change reaction into list format
  ## data process before using function
  rxnID <- unique(rxn$rxnID)
  rxn_annotation <- list()
  rxn_id <- list()
  for (i in 1:length(rxnID)) {
    rxn_id[[i]] <- which(rxn$rxnID %in% rxnID[i] == TRUE)
    rxn_annotation[[i]] <- rxn[rxn_id[[i]], ]
  }
  
  gpr0 <- list()
  for (i in 1:length(rxn_p0$rxnID)) {
    gpr0[[i]] <- rxn_p0[i, ]
  }
  
  ## using the main function
  txt6_new <- vector()
  for (j in 1:length(rxnID)) {
    txt6_new <- c(txt6_new, getRxnInformation(rxn_annotation[[j]], rxn_p = gpr0[[j]], R30=R3,R70=R7))
  }
  
  txt_new <- c(txt1, txt2, txt30, txt4, txt5, txt6_new, txt7) # connect each template
  
  # save the file
  # only with component coordinate
  writeLines(txt_final0, file("result/model_add_metabolite_protein.xml")) # save
  # adding the component connection
  writeLines(txt_new, file("result/model_with_protein-add rxn connection.xml")) # save
}
