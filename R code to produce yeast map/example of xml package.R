#example 1
library(XML)
doc = newXMLDoc()
document = newXMLNode("Document", doc = doc)
set = newXMLNode("Settings", parent = document)
elements = newXMLNode("Elements", parent = set)
newXMLNode("Canvas", parent = elements, attrs = c(Id = "1"))
newXMLNode("Canvas", parent = elements, attrs = c(Id = "2"))

objcol = newXMLNode("ObjectCollection", parent = document)
timeSeries1 = newXMLNode("Timeseries", parent = objcol)
timeSeries2 = newXMLNode("Timeseries", parent = objcol)

cat(saveXML(doc, 
            indent = TRUE, 
            prefix = "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"no\"?>\n"),
    file="example1.xml")



#example 1
library(XML)
doc = newXMLDoc()
document = newXMLNode("sbml", doc = doc)
set = newXMLNode("model", parent = document)
Notes = newXMLNode("notes", parent = set)
# main node
Annotation = newXMLNode("annotation", parent = set)
Celldesigner_extension = newXMLNode("celldesigner_extension",parent = Annotation)
Celldesigner_modelVersion = newXMLNode("celldesigner_modelVersion", parent = Celldesigner_extension)
Celldesigner_listOfSpeciesAliases = newXMLNode("celldesigner_listOfSpeciesAliases", parent = Celldesigner_extension)


#main node
ListOfUnitDefinitions = newXMLNode("listOfUnitDefinitions", parent = set)
UnitDefinition = newXMLNode("unitDefinition", parent = ListOfUnitDefinitions)
ListOfCompartments = newXMLNode("listOfCompartments", parent = set)
Compartments= newXMLNode("compartment", parent = ListOfCompartments )

ListOfSpecies = newXMLNode("listOfSpecies", parent = set)


cat(saveXML(doc, 
            indent = TRUE, 
            prefix = "<?xml version=\"1.0\" encoding=\"utf-8\">\n"),
    file="example1.xml")







#example 2
node = newXMLNode("A")
sapply(c("X", "Y", "Z", "X", "Y"),
       newXMLNode, parent = node)
cat(saveXML(node))
xmlAttrs(node)["src"] = "http://www.omegahat.org"



#example 3
library(XML)    
top = newXMLNode("a")
newXMLNode("b", attrs=c(x=1, y='abc'), parent=top)
newXMLNode("c", "With some text", parent=top)
top


#example 4
data <-
  read.csv(textConnection('"date","UYG.Open","UYG.High","UYG.Low","UYG.Close","UYG.Volume","UYG.Adjusted"
                          "2007-02-01",71.32,71.34,71.32,71.34,200,69.23
                          "2007-02-02",72.2,72.2,72.2,72.2,200,70.06
                          "2007-02-05",71.76,71.76,71.76,71.76,5100,69.63
                          "2007-02-06",72.85,72.85,72.85,72.85,3800,70.69
                          "2007-02-07",72.85,72.85,72.85,72.85,0,70.69'),
           as.is=TRUE)

library(XML)
xml <- xmlTree()
xml$addTag("document", close=FALSE)
for (i in 1:nrow(data)) {
  xml$addTag("row", close=FALSE)
  for (j in names(data)) {
    xml$addTag(j, data[i, j])
  }
  xml$closeTag()
}
xml$closeTag()

cat(saveXML(xml))
saveXML(xml, file="example3.xml")



#example 4
library(XML)

# XML STRING 
prefix.xml <- "<reports>
<report type='standard'>
<data> xxx </data>
<data> xxx </data>
<data> xxx </data>
</report>     
</reports>"

# DUMMY DATA FRAME
df <- data.frame("xxx","yyy")

# BUILD XML TREE
doc = xmlTreeParse(prefix.xml, useInternalNodes = T)     # PARSE STRING
root = xmlRoot(doc)                                      # FIND ROOT

reportNode = newXMLNode("report", parent=root)           # ADD TO ROOT
xmlAttrs(reportNode) = c(type = "enhanced")              # ADD ATTRIBUTE
pagesNode = newXMLNode("pages", parent=reportNode)       # ADD TO REPORT

for (i in 1:nrow(df)){
  pageNode = newXMLNode("page", parent=pagesNode)        # ADD PAGE FOR EACH RECORD
  for (j in 1:nrow(df)){
    newXMLNode("X.xxx.", df$X.xxx.[i], parent=pageNode)  # ADD COL/ROW VALUE
    newXMLNode("X.yyy.", df$X.yyy.[i], parent=pageNode)  # ADD COL/ROW VALUE
  }  
}

# VIEW XML
print(doc)

# SAVE XML TO FILE
saveXML(doc, file="example4.xml")





##example  5
data<- structure(list(page = c("Page One", "Page Two"), base64 = c("Very Long String thats a base64 character", "Very Long String thats a base64 character")), .Names = c("page", "base64"), row.names = 1:2, class = "data.frame")
names(data) <- c("title", "page")

library(XML)
xml <- xmlTree()
# names(xml)
xml$addTag("report", close=FALSE, attrs=c(type="enhanced"))
xml$addTag("pages", close=FALSE)
for (i in 1:nrow(data)) {
  xml$addTag("page", close=FALSE)
  for (j in names(data)) {
    xml$addTag(j, data[i, j])
  }
  xml$closeTag()
}
xml$closeTag()
xml$closeTag()
cat(saveXML(xml))
saveXML(xml, file="example5.xml")





