
## these code can be used the medium process of draft map after automatic layout
library(stringr)
library(tidyverse)
yeast_map <- readLines(file("model_test_protein-add rxn connection_using function_auto_layout.xml"))
anchor <- readLines(file("anchor.xml"))
index_linkAnchor <- which(str_detect(yeast_map,"celldesigner:linkAnchor")==FALSE)
yeast_map_new <- yeast_map[index_linkAnchor]

index_linkAnchor_arrow <- which(str_detect(yeast_map_new,"celldesigner:editPoints")==TRUE)
index_linkAnchor_arrow1 <- which(str_detect(yeast_map_new,"<celldesigner:lineDirection index=\"1\" value=\"unknown\"/>")==TRUE)
index_linkAnchor_arrow2 <- which(str_detect(yeast_map_new,"<celldesigner:lineDirection index=\"2\" value=\"unknown\"/>")==TRUE)
index_arrow0 <- c(index_linkAnchor_arrow,index_linkAnchor_arrow1,index_linkAnchor_arrow2)
yeast_map_new1 <- yeast_map_new[-index_arrow0]

writeLines(yeast_map_new1, file("model_test_protein-add rxn connection_using function_auto_layout_remove the linkAnchor.xml")) # save


##others: update the reaction reversibility

