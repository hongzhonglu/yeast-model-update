## find and replace____method 1
library(tidyverse)
library(xml2)
x <- read_xml("Acetate utilization.xml")
xml_find_all(x, "//annotation")
xml_structure(x, indent = 5)
write_xml(x, "Acetate utilization 2.xml", options = "format")



## find and replace______method 2
txt <- readLines(file("Acetate utilization.xml"))
##replace the metabolites
m <- which(str_detect(txt, "<listOfSpecies>") == TRUE)
n <- which(str_detect(txt, "</listOfSpecies>") == TRUE)
for (i in m:n){
  txt[i] <- str_replace_all(txt[i],"name=\"Acetate\"","name=\"Acetate_new\"")
  txt[i] <- str_replace_all(txt[i],"AcetateSynonym","Acetate_newSynonym")
  txt[i] <- str_replace_all(txt[i],">Acetate<",">Acetate_new<")
  
  }

writeLines(txt, file("test.xml"))



## STEP 1
## if we only add a new metabolites but without a reaction
##<celldesigner:speciesAlias id="sa7" species="s12">
  ##<celldesigner:activity>inactive#</celldesigner:activity>
  #<celldesigner:bounds x="24.0" y="25.5" w="70.0" h="25.0"/>
  #<celldesigner:font size="12"/>
  #<celldesigner:view state="usual"/>
  #<celldesigner:usualView>
  #<celldesigner:innerPosition x="0.0" y="0.0"/>
  #<celldesigner:boxSize width="70.0" height="25.0"/>
  #<celldesigner:singleLine width="1.0"/>
  #<celldesigner:paint color="ffccff66" scheme="Color"/>
  #</celldesigner:usualView>
  #<celldesigner:briefView>
  #<celldesigner:innerPosition x="0.0" y="0.0"/>
  #<celldesigner:boxSize width="80.0" height="60.0"/>
  #<celldesigner:singleLine width="0.0"/>
  #<celldesigner:paint color="3fff0000" scheme="Color"/>
  #</celldesigner:briefView>
  #<celldesigner:info state="empty" angle="-1.5707963267948966"/>
  #</celldesigner:speciesAlias>


  #<species metaid="s12" id="s12" name="AMP" compartment="default" initialAmount="0">
  #<annotation>
  #<celldesigner:extension>
  #<celldesigner:positionToCompartment>inside#</celldesigner:positionToCompartment>
  #<celldesigner:speciesIdentity>
  #<celldesigner:class>SIMPLE_MOLECULE#</celldesigner:class>
  #<celldesigner:name>AMP#</celldesigner:name>
  #</celldesigner:speciesIdentity>
  #</celldesigner:extension>
  #</annotation>
  #</species>


##STEP 2
## we add a reaction based on the metabolites
#<reaction metaid="re3" id="re3" reversible="false">
  #<annotation>
  #<celldesigner:extension>
  #<celldesigner:reactionType>STATE_TRANSITION#</celldesigner:reactionType>
  #<celldesigner:baseReactants>
  #<celldesigner:baseReactant species="s9" alias="a10"/>
  #</celldesigner:baseReactants>
  #<celldesigner:baseProducts>
  #<celldesigner:baseProduct species="s11" alias="sa6">
  #<celldesigner:linkAnchor position="S"/>
  #</celldesigner:baseProduct>
  #</celldesigner:baseProducts>
  #<celldesigner:connectScheme connectPolicy="direct" rectangleIndex="0">
  #<celldesigner:listOfLineDirection>
  #<celldesigner:lineDirection index="0" value="unknown"/>
  #</celldesigner:listOfLineDirection>
  #</celldesigner:connectScheme>
  #<celldesigner:line width="1.0" color="ff000000"/>
  #</celldesigner:extension>
  #</annotation>
  #<listOfReactants>
  #<speciesReference metaid="CDMT00001" species="s9">
  #<annotation>
  #<celldesigner:extension>
  #<celldesigner:alias>a10#</celldesigner:alias>
  #</celldesigner:extension>
  #</annotation>
  #</speciesReference>
  #</listOfReactants>
  #<listOfProducts>
  #<speciesReference metaid="CDMT00007" species="s11">
  #<annotation>
  #<celldesigner:extension>
  #<celldesigner:alias>sa6#</celldesigner:alias>
  #</celldesigner:extension>
  #</annotation>
  #</speciesReference>
  #</listOfProducts>
  #</reaction>



## step 3
## if we want to add a new metabolite as well as  a new reaction
## we should update
##  model > annotation >  cellDesigner:listOfSpeciesAliases
##  listOfSpecies > species s11
##  listOfReactions > reaction
##  listOfReactions > listOfProduct


##  id of metabolites
##  <celldesigner:speciesAlias id="sa6" species="s11">
##  metaid="s11" id="s11" name="AMP" 
##  metaid="CDMT00007" species="s11">


## if we only add a new reaction but without new metabolites
##<reaction metaid="re4" id="re4" reversible="false">
##<annotation>
##<celldesigner:extension>
#<celldesigner:reactionType>STATE_TRANSITION#</celldesigner:reactionType>
#<celldesigner:baseReactants>
#<celldesigner:baseReactant species="s2" alias="a11">
#<celldesigner:linkAnchor position="WNW"/>
#</celldesigner:baseReactant>
#</celldesigner:baseReactants>
#<celldesigner:baseProducts>
#<celldesigner:baseProduct species="s11" alias="sa6">
#<celldesigner:linkAnchor position="E"/>
#</celldesigner:baseProduct>
#</celldesigner:baseProducts>
#<celldesigner:connectScheme connectPolicy="direct" rectangleIndex="0">
#<celldesigner:listOfLineDirection>
#<celldesigner:lineDirection index="0" value="unknown"/>
#</celldesigner:listOfLineDirection>
#</celldesigner:connectScheme>
#<celldesigner:line width="1.0" color="ff000000"/>
#</celldesigner:extension>
#</annotation>
#<listOfReactants>
#<speciesReference metaid="CDMT00008" species="s2">
#<annotation>
#<celldesigner:extension>
#<celldesigner:alias>a11#</celldesigner:alias>
#</celldesigner:extension>
#</annotation>
#</speciesReference>
#</listOfReactants>
#<listOfProducts>
#<speciesReference metaid="CDMT00009" species="s11">
#<annotation>
#<celldesigner:extension>
#<celldesigner:alias>sa6#</celldesigner:alias>
#</celldesigner:extension>
#</annotation>
#</speciesReference>
#</listOfProducts>
#</reaction>
#

##step 4
#if we want to add a reactions with more than two metabolites



