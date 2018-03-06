library(stringr)

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
  txt [16] <-  paste("<celldesigner:modelDisplay sizeX=\"",x,"\" sizeY=\"",y,"\"/>",sep="")                                                                                                                                                                                                  
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


### test 1  add metabolites using function from the imported data
library(readxl)
met <- read_excel("Metbolites.xlsx", sheet = "Sheet1")
id <- met$id
species <- met$species
name <- met$name
x <- as.character(met$x)
y <- as.character(met$y)

txt2 <- vector()
txt4 <- vector()
for (i in seq(length(id))){
  txt2 <- c(getPosition(id01=id[i],species01=species[i],x01=x[i],y01=y[i]), txt2)
}

species_unique <- unique(met$species)
name_unique <- unique(met$name)

for (i in seq(length(species_unique))){
  txt4 <- c(getAnnotation(species01 = species_unique[i],name01 = name_unique[i]), txt4)
}

txt_final <- c(txt1,txt2,txt3,txt4,txt5,txt7)
writeLines(txt_final, file("template_add_metabolite_test.xml")) # save




### Application 1 import reactions from transport subsystem
library(readxl)
met_transport0<- read_excel("metabolit_from_transport_yeastGEM.xlsx")

##function to give the postion information of each metabolite

zone_rxn <- read_excel("metabolit_from_transport_yeastGEM.xlsx", 
           sheet = "Sheet2")


#x position
zone_rxn$rxnID[1]
zone_rxn$x[1]
zone_rxn$y[1]

r_position <- list()
p_position <- list()
reactant_len <- vector()
product_len <- vector()
y_distance1 <- vector()
y_distance2 <- vector()
for (i in 1:100){
r_position[[i]] <- which(met_transport0$rxnID==zone_rxn$rxnID[i] & met_transport0$type =="reactant")
p_position[[i]] <- which(met_transport0$rxnID==zone_rxn$rxnID[i] & met_transport0$type =="product")
met_transport0$x[r_position[[i]]] <- zone_rxn$x[i]-250
met_transport0$x[p_position[[i]]] <- zone_rxn$x[i]


#y position
reactant_len[i] <- length(r_position[[i]])
product_len[i] <- length(p_position[[i]])
y_distance1[i] <- 200/reactant_len[i]
y_distance2[i] <- 200/product_len[i]

met_transport0$y[r_position[[i]]] <- seq(from = zone_rxn$y[i]-200, to = zone_rxn$y[i], by = y_distance1[i])[1:reactant_len[i]] 
met_transport0$y[p_position[[i]]] <- seq(from = zone_rxn$y[i]-200, to = zone_rxn$y[i], by = y_distance2[i])[1:product_len[i]]

}


met_transport0 <- met_transport0[1:336,]

id <- met_transport0$id
species <- met_transport0$species
name <- met_transport0$name
x <- as.character(met_transport0$x)
y <- as.character(met_transport0$y)

txt1 <-getNewSize(x=4000,y=2000)
txt2 <- vector()
txt4 <- vector()
for (i in seq(length(id))){
  txt2 <- c(getPosition(id01=id[i],species01=species[i],x01=x[i],y01=y[i]), txt2)
}

species_unique <- unique(met_transport0$species)
name_unique <- unique(met_transport0$name)

for (i in seq(length(species_unique))){
  txt4 <- c(getAnnotation(species01 = species_unique[i],name01 = name_unique[i]), txt4)
}

txt_final <- c(txt1,txt2,txt3,txt4,txt5,txt7)
writeLines(txt_final, file("model_add_metabolite_test.xml")) # save




### Application2 import reactions from core carbon metabolism EMP, PPP and TCA cycle
library(readxl)
met_transport0<- read_excel("metabolit_from_core_carbon_yeastGEM.xlsx")

##function to give the postion information of each metabolite

zone_rxn <- read_excel("metabolit_from_core_carbon_yeastGEM.xlsx", 
                       sheet = "Sheet2")


#x position
zone_rxn$rxnID[1]
zone_rxn$x[1]
zone_rxn$y[1]

for (i in 1:length(zone_rxn$rxnID)){
  met_transport0$x[i] <- zone_rxn$x[i]
  met_transport0$y[i] <- zone_rxn$y[i]
}


met_transport0 <- met_transport0[1:length(met_transport0$id),]

id <- met_transport0$id
species <- met_transport0$species
name <- met_transport0$name
x <- as.character(met_transport0$x)
y <- as.character(met_transport0$y)

txt1 <-getNewSize(x=2200,y=1200)
txt2 <- vector()
txt4 <- vector()
for (i in seq(length(id))){
  txt2 <- c(getPosition(id01=id[i],species01=species[i],x01=x[i],y01=y[i]), txt2)
}

species_unique <- unique(met_transport0$species)
name_unique <- unique(met_transport0$name)

for (i in seq(length(species_unique))){
  txt4 <- c(getAnnotation(species01 = species_unique[i],name01 = name_unique[i]), txt4)
}

txt_final <- c(txt1,txt2,txt3,txt4,txt5,txt7)
writeLines(txt_final, file("model_add_metabolite_test1.xml")) # save






