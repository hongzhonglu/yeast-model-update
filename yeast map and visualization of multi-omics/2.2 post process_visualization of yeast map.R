model_visualization <- readLines(file("data/model_visualization.xml"))
#gpr_position_in_map <- gpr

# extabolish the position of gene or protein

"<celldesigner:speciesAlias id=\"sa8\" species=\"s8\">"
gene_list <- "<celldesigner:speciesAlias id=\"sa8\" species=\"s8\">"
color_position_g <- which(str_detect(model_visualization, gene_list)==TRUE)+9 #fill color

#color_position_g <- which(str_detect(model_visualization, gene_list)==TRUE)+15 #line color
model_visualization[color_position_g] <- "<celldesigner:paint color=\"ff0eb10e\" scheme=\"Color\"/>"   #green

model_visualization[color_position_g] <- "<celldesigner:paint color=\"ffff0000\" scheme=\"Color\"/>"   #green


# estabolish the postion of reaction

"<reaction metaid=\"r_0200\" id=\"r_0200\" reversible=\"false\">"
size_reaction <- which(str_detect(model_visualization,"<celldesigner:line width=\"1.0\" color=\"ff000000\"/>")==TRUE)
model_visualization[size_reaction[1]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[4]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[7]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[10]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[13]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[16]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"





writeLines(model_visualization, file("result/model_visualization_changed_color.xml")) # save

