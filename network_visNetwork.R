##### functions --------------------

#' helper function to install and load dependencies
#' 
ortSuite_depend <- function(){
  # not installed 
  not_i <- setdiff(
    c("visNetwork","reshape2","RColorBrewer"), 
    rownames(installed.packages())
  )
  # install packages
  if(length(not_i)>0){
    install.packages(not_i)  
  }
  # library
  library(visNetwork)
  library(reshape2)
  library(RColorBrewer)
}

#' OrtSuite Network visualization 
#' 
#' @param species_reaction Species (columns) x reactions (rows) table
#' @param reaction_reaction Reaction reaction connection list (2 columns)
#' @param color Color nodes and species. Possible values are reaction, species and type. Reaction: each 
#' reaction node and associated species node get different color. Species: each species node get a 
#' different color. Type: species and reactions get different colors.
ortSuite_net <- function(species_reaction, reaction_reaction = NULL, color = 'reaction'){
  
  ## filter reactions
  # if reaction_reaction list given, filter by this
  if(!is.null(reaction_reaction)){
    reaction_u <- unique(c(reaction_reaction[,1],reaction_reaction[,2]))
  # if reaction_reaction not given, take all reactions from species_reaction table
  } else {
    reaction_u <- unique(rownames(species_reaction))
  }
  
  # filter for valid reactions
  species_reaction_val <- species_reaction[is.element(rownames(species_reaction),reaction_u),]
  
  # create edge list (species - reactions)
  species_reaction_list <- melt(cbind(rownames(species_reaction_val),species_reaction_val),id=1)
  # make character (from factor)
  for(i in 1:2){
    species_reaction_list[,i] <- as.character(species_reaction_list[,i])
  }
  species_reaction_list <- species_reaction_list[species_reaction_list[,3]==1,1:2]
  colnames(species_reaction_list) <- c('from', 'to')
  # store unique  species (for later use)
  spec_u <- unique(species_reaction_list[,2])
  # for each reaction, create a virtual copy of the species
  reac_u <- unique(species_reaction_list$from)
  for(i in 1:length(reac_u)){
    # make species unique 
    species_reaction_list[species_reaction_list[,1]==reac_u[i],2] <- paste0(species_reaction_list[species_reaction_list[,1]==reac_u[i],2],'_',i)
  }
  
  # add reaction - reaction edges
  species_reaction_list <- rbind(
    species_reaction_list,
    reaction_reaction
  )
  
  ### build node list
  # ids
  node_ids <- unique(c(species_reaction_list[,1],species_reaction_list[,2]))
  # shape
  node_shape <- rep('square',length(node_ids)) # set all to square
  node_shape[grep("_[0-9]{1,3}$",node_ids)] <- 'dot' # species to dot (small circle)
  ### color
  if(color == 'type'){
    # color by type (reaction and species)
    node_color <- rep('darkred',length(node_ids))
    node_color[grep("_[0-9]{1,3}$",node_ids)] <- 'dodgerblue'
  } else if(color == 'species') {
    # color by species
    color_set <- colorRampPalette(brewer.pal(8, "Dark2"))(length(spec_u))
    node_color <- rep('darkred',length(node_ids))
    for(i in 1:length(spec_u)){
      node_color[is.element(gsub("_[0-9]{1,3}$","",node_ids),spec_u[i])] <- color_set[i]
    }
  } else if(color == 'reaction') {
    # color by reaction
    color_set <- colorRampPalette(brewer.pal(8, "Dark2"))(length(reac_u))
    node_color <- rep('black',length(node_ids))
    for(i in 1:length(reac_u)){
      node_color[grep(paste0("_",i,"$"),node_ids)] <- color_set[i] # color assiciated species
      node_color[grep(reac_u[i],node_ids)] <- color_set[i] # color reaction
    }
  } else {
    stop('Wrong color parameter given')
  }
  # tooltips
  node_title <- paste0('<p>Reaction <b>',node_ids,'</b><br><a href="https://www.genome.jp/dbget-bin/www_bget?rn:',node_ids,'">KEGG link</a></p>')
  node_title[grep("_[0-9]{1,3}$",node_ids)] <- paste0('<p>Species <b>',gsub("_[0-9]{1,3}$","",node_ids[grep("_[0-9]{1,3}$",node_ids)]),'</b></p>')
  # node size
  node_size <- rep(1,length(node_ids))
  node_size[grep("_[0-9]{1,3}$",node_ids)] <- 1
  
  # create node data frame
  nodes <- data.frame(
    # node id
    id = node_ids,
    # label 
    label = paste(gsub("_[0-9]{1,3}$","",node_ids)),
    # tooltip (HTML)
    title = node_title,
    # size
    value = node_size,
    # shape
    shape = node_shape,
    # color
    color = node_color,
    # shadow 
    shadow = F,
    # size
    font.size=30
  )
  
  # edges
  edges <-  species_reaction_list
  edges$color <- 'black'
  
  # build legend
  num_col <- 1
  if(color == 'type'){
    node_legend <- data.frame(
      color = c('darkred','dodgerblue'),
      shape = c('square','dot'),
      label = c('Reaction','Species')
    )
  } else if(color == 'reaction'){
    node_legend <- data.frame(
      color = node_color[is.element(node_ids,reaction_u)],
      shape = 'square',
      label = node_ids[is.element(node_ids,reaction_u)]
    )
  } else if(color == 'species'){
    node_legend <- data.frame(
      color = color_set,
      shape = 'dot',
      label = spec_u
    )
  }
  # increase number of legend columns
  if(nrow(node_legend)>10){
    num_col <- 2
  } 
  
  # create network
  visNetwork(nodes, edges = edges, width = '100%', height = '800px') %>%
    visLegend(addNodes = node_legend, useGroups = F, ncol = num_col) %>%
    visOptions(selectedBy = "label") 
}

##### Create network --------------------

### get parameters
get_param <- commandArgs(trailingOnly = TRUE)

# check input
if(length(get_param)<1){
  stop('No input files')
}

### load and parse files
# load species reaction table 
s_r <- read.table(
  file = get_param[1], 
  header = T, sep = ",",
  stringsAsFactors = F,
  row.names = 1
)

# load reaction (to reaction) table
if(!is.na(get_param[2])){
  r_r <- read.table(
    file = get_param[2], 
    header = F, sep = "\t",
    stringsAsFactors = F
  )
  colnames(r_r) <- c('from','to')
} else {
  r_r <- NULL
}
### install and load dependencies
ortSuite_depend()
  
### create and save network
net <- ortSuite_net(s_r, r_r)
visSave(net,'network_ortSuite.html')

