library(stringr)
library(tidyverse)
## connect the different template
txt1 <- readLines(file("template1")) # model information and size
#txt2 <- readLines(file("template2")) # repeated part: position information of each metabolites
txt3 <- readLines(file("template3")) # other information of model besides metabolits and reactions
#txt4 <- readLines(file("template4")) # repeated part: annotation of each metabolites; should be unique id
txt5 <- readLines(file("template5")) # sign of xml format
#txt6 <- readLines(file("template6")) # repeated part: information of each reactions
txt7 <- readLines(file("template7")) # end sign of xml format


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

readLines(file("txt4_gene annotation template.xml"))


# define annotation information for each gene
getGeneAnnotation <- function(species01,name01, geneID01){
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











###application 1 import metabolite from core carbon metabolism

### Application2 import reactions from core carbon metabolism EMP, PPP and TCA cycle
library(readxl)
#example:met_transport0<- read_excel("metabolit_from_core_carbon_yeastGEM.xlsx")

txt1 <-getNewSize(x=2500,y=3600)
txt2 <- vector()
txt3 <- vector()
txt4 <- vector()


# txt2
id <- met_annotation$id
species <- met_annotation$species
name <- met_annotation$name
x <- as.character(met_annotation$x)
y <- as.character(met_annotation$y)

for (i in seq(length(id))){
  txt2 <- c(getPosition(id01=id[i],species01=species[i],x01=x[i],y01=y[i]), txt2)
}


# txt3
## add protein define
## add gene define
#example: protein <- read_excel("Metbolites.xlsx", sheet = "protein")
#example: gene <- read_excel("Metbolites.xlsx", sheet = "gene")

## input protein and gene information
protein0 <- filter(met_annotation, type=="PROTEIN")
protein0$type <- "GENERIC"
protein0$protein_id <- paste("pr",1:length(protein0$metID), sep="")

gene0 <- filter(met_annotation, type=="GENE")
gene0$gene_id <- paste("gn",1:length(gene0$metID), sep="")


getProteinGene <- function (id,type,name){
  txt0 <- vector()
  if (type =="GENERIC"){
    txt0[1] <- paste("<celldesigner:protein id=\"",id,"\" name=\"",name,"\" type=\"GENERIC\"/>", sep="")
  }
  if(type =="GENE"){
    txt0[1] <- paste("<celldesigner:gene id=\"",id,"\" name=\"",name,"\" type=\"GENE\"/>", sep="")     
  }
  
  return(txt0)
}


txt3_1 <- readLines(file("txt3_1.xml"))

#txt3_2 protein module
protein_id <- protein0$protein_id
type_p <- protein0$type
name_p <- protein0$name


txt3_2 <- vector()
for(i in seq(length(protein_id))){
  txt3_2 <- c(txt3_2, getProteinGene(id=protein_id[i], type = type_p[i], name = name_p[i]))
}


txt3_3 <- readLines(file("txt3_3.xml"))

#txt3_4 gene module
gene_id <- gene0$gene_id
type_g <- gene0$type
name_g <- gene0$name
txt3_4 <- vector()
for (i in seq(length(gene_id))){
  txt3_4 <- c(txt3_4, getProteinGene(id=gene_id[i], type = type_g[i], name=name_g[i]))
}


txt3_5 <- readLines(file("txt3_5.xml"))
txt30 <- c(txt3_1,txt3_2,txt3_3,txt3_4,txt3_5)


#txt4
##txt4: give the annotation for each unique metabolite
#txt4_1 metabolite module
#example: met <- read_excel("Metbolites.xlsx", sheet = "metabolite")
#example: met_without_protein <- filter(met,type == "SIMPLE_MOLECULE")

met_without_protein0 <- filter(met_annotation, type == "SIMPLE_MOLECULE")
species_unique <-unique(met_without_protein0$species)
name_unique <- unique(met_without_protein0$name)

txt4_1 <- vector()
for (i in seq(length(species_unique))){
  txt4_1 <- c(txt4_1, getAnnotation(species01 = species_unique[i], name01 = name_unique[i]))
}


#txt4_2 protein module
#example: protein <- read_excel("Metbolites.xlsx", sheet = "protein")

protein0 <- filter(met_annotation,type == "PROTEIN")
protein0$protein_id <- paste("pr",1:length(protein0$metID), sep="")
protein_specie <- protein0$species
protein_id <- protein0$protein_id
type_p <- protein0$type
name_p <- protein0$name

txt4_2 <- vector()
for(i in seq(length(protein_id))){
  txt4_2 <- c(txt4_2, getProAnnotation(species01 = protein_specie[i],name01 = name_p[i], proteinID01=protein_id[i]))
}


#txt4_3 gene module
gene_id <- gene0$gene_id
gene_specie <- gene0$species
type_g <- gene0$type
name_g <- gene0$name
txt4_3 <- vector()
for (i in seq(length(gene_id))){
  txt4_3 <- c(txt4_3, getGeneAnnotation(species01=gene_specie[i],name01=name_g[i], geneID01=gene_id[i]))
}
txt4 <- c(txt4_1,txt4_2,txt4_3)

txt_final0 <- c(txt1,txt2,txt30,txt4,txt5,txt7)
writeLines(txt_final0, file("model_add_metabolite_protein.xml")) # save




















### test 1  add metabolites using function from the imported data
library(readxl)
met <- read_excel("Metbolites.xlsx", sheet = "metabolite")

#txt1 : define the size of graph
txt1 <- getNewSize(x=1200, y=1200)

#txt2 : define the postion information for the metabolites, proteins and genes
id <- met$id
species <- met$species
name <- met$name
met_type <- met$type

x <- as.character(met$x)
y <- as.character(met$y)

txt2 <- vector()

for (i in seq(length(id))){
  txt2 <- c(getPosition(id01=id[i],species01=species[i],x01=x[i],y01=y[i]), txt2)
}


## txt3 Replace in the above module
## add protein define
## add gene define
protein <- read_excel("Metbolites.xlsx", sheet = "protein")
gene <- read_excel("Metbolites.xlsx", sheet = "gene")
getProteinGene <- function (id,type,name){
  txt0 <- vector()
  if (type =="GENERIC"){
    txt0[1] <- paste("<celldesigner:protein id=\"",id,"\" name=\"",name,"\" type=\"GENERIC\"/>", sep="")
  }
  if(type =="GENE"){
    txt0[1] <- paste("<celldesigner:gene id=\"",id,"\" name=\"",name,"\" type=\"GENE\"/>", sep="")     
  }
  
  return(txt0)
}


txt3_1 <- readLines(file("txt3_1.xml"))

#txt3_2 protein module
protein_id <- protein$protein_id
type_p <- protein$type
name_p <- protein$name
txt3_2 <- vector()
for(i in seq(length(protein_id))){
  txt3_2 <- c(txt3_2, getProteinGene(id=protein_id[i], type = type_p[i], name = name_p[i]))
}


txt3_3 <- readLines(file("txt3_3.xml"))

#txt3_4 gene module
gene_id <- gene$gene_id
type_g <- gene$type
name_g <- gene$name
txt3_4 <- vector()
for (i in seq(length(gene_id))){
  txt3_4 <- c(txt3_4, getProteinGene(id=gene_id[i], type = type_g[i], name=name_g[i]))
}

txt3_5 <- readLines(file("txt3_5.xml"))

txt30 <- c(txt3_1,txt3_2,txt3_3,txt3_4,txt3_5)

##txt4: give the annotation for each unique metabolite
#txt4_1 metabolite module
met_without_protein <- filter(met,type == "SIMPLE_MOLECULE")
unique <- unique(met_without_protein$species) ####should be noted for this function
index0 <- which(met_without_protein$species %in% unique ==TRUE)
met0 <- met_without_protein[index0,]
species_unique <-met0$species
name_unique <- met0$name
class_type <- met0$type
txt4_1 <- vector()
for (i in seq(length(species_unique))){
  txt4_1 <- c(txt4_1, getAnnotation(species01 = species_unique[i], name01 = name_unique[i]))
}

#txt4_2 protein module
protein <- read_excel("Metbolites.xlsx", sheet = "protein")
protein_specie <- protein$species
protein_id <- protein$protein_id
type_p <- protein$type
name_p <- protein$name
txt4_2 <- vector()
for(i in seq(length(protein_id))){
  txt4_2 <- c(txt4_2, getProAnnotation(species01 = protein_specie[i],name01 = name_p[i], proteinID01=protein_id[i]))
}

#txt4_3 gene module
gene_id <- gene$gene_id
gene_specie <- gene$species
type_g <- gene$type
name_g <- gene$name
txt4_3 <- vector()
for (i in seq(length(gene_id))){
  txt4_3 <- c(txt4_3, getGeneAnnotation(species01=gene_specie[i],name01=name_g[i], geneID01=gene_id[i]))
}
txt4 <- c(txt4_1,txt4_2,txt4_3)

txt_final0 <- c(txt1,txt2,txt30,txt4,txt5,txt7)
writeLines(txt_final0, file("template_add_metabolite_protein.xml")) # save









