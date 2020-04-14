

args <- commandArgs(trailingOnly = TRUE)  

library("readxl")
library(stringr)
library(dplyr)
library(gsubfn)
library(mgsub)
library(stringr)
library(reshape)
library(reticulate)


#output_ortan <- read.csv("~/OrtSuite/examples/OrtAn_Results/Results/Species_Annotation.csv",header=T,sep=";",row.names = 1) # Read in output file from OrtAn
#gpr_file <- read_xlsx("~/OrtSuite/examples/test_database/gpr.xlsx") # read in gpr file generated from .sh script
#user_input <- read.csv("~/OrtSuite/examples/OrtAn_Results/Results/test_user_input.csv",header = T) # read in user input file for constraints


  output_ortan <- read.csv(args[2],header=T,sep=";",row.names = 1) # Read in output file from OrtAn
  gpr_file <- read_xlsx(args[1]) # read in gpr file generated from .sh script
  user_input <- read.csv(args[3],header = T) # read in user input file for constraints

## Create files needed for subsequent tasks

temp_reaction_table <- unique(gpr_file$Reaction) # Get unique reaction ids

final_reaction_table <- matrix(NA,nrow=length(temp_reaction_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_reaction_table) <- temp_reaction_table
colnames(final_reaction_table) <- "Definition"

# paste all possible gpr rules for each unique reaction id

for (reaction in 1:length(temp_reaction_table)){
  x <- gpr_file[gpr_file$Reaction %in% temp_reaction_table[reaction],]
  
  final_reaction_table[reaction,] <- paste(x$Definition, collapse = ";", sep=",")

}


#Check if multiple rows for the same reaction id exist. If so check if reaction ids exist are occuring in multiple lines and keep only the second column.
reaction_table <- matrix(NA,nrow=nrow(final_reaction_table),ncol=1)

for (line in 1:nrow(final_reaction_table)){
  tt <- strsplit(final_reaction_table[line,1],";")
  if(grepl(tt[[1]][1],tt[[1]][2]))
     reaction_table[line,1]=tt[[1]][2]
  else
    reaction_table[line,1]=final_reaction_table[line,1]
  
  
  row.names(reaction_table) <- rownames(final_reaction_table)
  colnames(reaction_table) <- "Definition"
  
  
}

#Replace all "and" by [reaction a, reaction b].
#Replace all or by "," without the square brackets.

for(line in 1:nrow(reaction_table)){
  if(str_detect(reaction_table[line,1],"and"))
    reaction_table[line,1] <- paste("[",gsub("and",",",reaction_table[line,1]),"]")
  if(str_detect(reaction_table[line,1],"or"))
    reaction_table[line,1] <- gsub("or",",",reaction_table[line,1])
  else
    reaction_table[line,1] <- reaction_table[line,1]
    }
# Output file ---> GPR rules for each reaction
'
       Definition            
R02601 "K18364 ,  K02554"    
R00750 "[ K18365 ,  K01666 ]"
R00813 "K05783"              
R02604 "K01617 ,  K10216"    
R00816 "K00446 ,  K07104" 

'


############ Obtain kos associated with each EC number ################

# Data obtained from gpr_file

temp_ec_table <- unique(gpr_file$Protein) # Get unique EC numbers

final_ec_table <- matrix(NA,nrow=length(temp_ec_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_ec_table) <- temp_ec_table
colnames(final_ec_table) <- "kos"

# paste all possible gpr rules for each unique EC number

for (ec in 1:length(temp_ec_table)){
  x <- gpr_file[gpr_file$Protein %in% temp_ec_table[ec],]
  
  final_ec_table[ec,] <- paste(x$Definition, collapse = ";", sep=",")
}


#Check if multiple rows for the same reaction id exist. If so check if reaction ids exist are occuring in multiple lines and keep only the second column.
ec_table <- matrix(NA,nrow=nrow(final_ec_table),ncol=1)

for (line in 1:nrow(final_ec_table)){
  tt <- strsplit(final_ec_table[line,1],";")
  if(grepl(tt[[1]][1],tt[[1]][2]))
    ec_table[line,1]=tt[[1]][2]
  else
    ec_table[line,1]=final_ec_table[line,1]
  
  
  row.names(ec_table) <- rownames(final_ec_table)
  colnames(ec_table) <- "GPR"
  
  
}

#Replace all "and" by [reaction a, reaction b].
#Replace all or by "," without the square brackets.

for(line in 1:nrow(ec_table)){
  if(str_detect(ec_table[line,1],"and"))
    ec_table[line,1] <- paste("[",gsub("and",",",ec_table[line,1]),"]")
  if(str_detect(ec_table[line,1],"or"))
    ec_table[line,1] <- gsub("or",",",ec_table[line,1])
  else
    ec_table[line,1] <- ec_table[line,1]
}


#Output file  ---> GP rules for each enzyme (will be used when transforming to .json)
'
          GPR                   
4.2.1.80  "K18364 ,  K02554"    
4.1.3.39  "[ K18365 ,  K01666 ]"
1.3.1.25  "K05783"              
3.7.1.9   "K01617 ,  K10216"    
1.13.11.2 "K00446 ,  K07104" 

'



############# Obtain reaction associated with EC numbers ##############

# Data obtained from gpr_file

temp_reac_ec_table <- unique(gpr_file$Reaction) # Get unique reaction ids

final_reac_ec_table <- matrix(NA,nrow=length(temp_reac_ec_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_reac_ec_table) <- temp_reac_ec_table
colnames(final_reac_ec_table) <- "Definition"

# paste all possible EC numbers for each unique reaction id

for (reaction in 1:length(temp_reac_ec_table)){
  x <- gpr_file[gpr_file$Reaction %in% temp_reac_ec_table[reaction],]
  
  final_reac_ec_table[reaction,] <- paste(x$Protein, collapse = ";", sep=",")
  
}


#Check if multiple rows for the same reaction id exist. If so check if reaction ids exist are occuring in multiple lines and keep only the second column.
reac_ec_table <- matrix(NA,nrow=nrow(final_reac_ec_table),ncol=1)

for (line in 1:nrow(final_reac_ec_table)){
  tt <- strsplit(final_reac_ec_table[line,1],";")
  if(grepl(tt[[1]][1],tt[[1]][2]))
    reac_ec_table[line,1]=tt[[1]][2]
  else
    reac_ec_table[line,1]=final_reac_ec_table[line,1]
  
  
  row.names(reac_ec_table) <- rownames(final_reac_ec_table)
  colnames(reac_ec_table) <- "EC_number"
  
  
}


#Output file ---> EC numbers associated with each reaction
'
      EC_number  
R02601 "4.2.1.80" 
R00750 "4.1.3.39" 
R00813 "1.3.1.25" 
R02604 "3.7.1.9"  
R00816 "1.13.11.2" 
'




####################### Create table mapping reactions to species based on gpr rules  ################################################


#Create new matrix where reactions will be stored mapped to species
reaction_species <- matrix(NA,nrow = nrow(reaction_table),ncol=ncol(output_ortan))

rownames(reaction_species) <- rownames(reaction_table)
colnames(reaction_species) <- colnames(output_ortan)

for(i in 1:nrow(reaction_species)){
  kos <- reaction_table[i,1] #get the GPr rules from the reaction table for line i
  for(j in 1:ncol(output_ortan)){ # iterate through each species in output_ortan table
    species_results <- as.matrix(as.data.frame(output_ortan[,j],row.names = rownames(output_ortan))) # get only a single column j
    if(grepl("\\[",kos)==T){ # check if there are "[" in the cell. Square brackets indicate that the set ok KOs must be all present in order to perform the reaction
      kos <- gsub("\\[","",kos) #replace [
      kos <- gsub("\\]","",kos) #replace ]
      kos <- gsub(" ","",kos) # remove white spaces
      x <- strsplit(kos,",") # split the string into n parts
      if(sum(species_results[which(rownames(species_results) %in% x[[1]])]) == length(x[[1]])){ #the sum must be equal to the number of KOs in gpr rule
        reaction_species[i,j] <- 1 # if all KOs are present then assign 1
      }else{
        reaction_species[i,j] <- 0
      }
    }else{
      kos <- gsub(" ","",kos) # remove empty spaces
      x <- strsplit(kos,",") # split string into n parts
      if(sum(species_results[which(rownames(species_results) %in% x[[1]])]) >= 1){ # at least one KO must have a value of one in species results table
        reaction_species[i,j] <- 1
      }else{
        reaction_species[i,j] <- 0
      }
    }
  }
}

#  Output file ---> A table mapping reactions to species already taking into account GR rules
'
       adv.faa ath.faa aza.faa azd.faa azi.faa test test_complete_1
R02601       1       0       0       1       1    1               1
R00750       1       0       0       1       1    1               1
R00813       0       0       0       0       0    1               1
R02604       1       0       0       1       0    0               1
R00816       1       0       0       1       0    0               1
R00001       0       1       0       1       1    0               0
'

write.csv(reaction_species,file="Reactions_mapped_to_species.txt")

#### Table with user defined constraints ########



extract_total_reactions <- strsplit(as.character(user_input$Reactions),",") # returns a list of 2 (in this case) sets
extract_subset_single_org <- str_match_all(user_input$Single_org, "(?<=\\().+?(?=\\))") # returns a list of 2 sets each containing the subsets of reations that need to be performed by a single species
extract_transporters <- str_match_all(user_input$Transporter, "(?<=\\().+?(?=\\))") # returns a list of 2 sets each containing the subsets of transport reactions that are needed for the subsets of single_org
names(extract_transporters) <- user_input$Pathway # Assign pathway names to transporter sets

#### Script to obtain all species individually capable of performing complete pathways ##########

complete_pathways <- list() # output list for species with complete pathways


list_names <- as.character(user_input$Pathway) # save pathway names from user_input

for(i in 1:length(extract_total_reactions)){ # for each pathway reaction set
  if(length(extract_transporters[[i]])==0){
    
  temp_species <- matrix(ncol=1) #temp file to store species that have complete pathway
  
  temp_total_reactions <- as.matrix(extract_total_reactions[[i]]) # temp file to store reactions for complete pathway
  
  for(j in 1:ncol(reaction_species)){ # iterate through each species in reaction_species table
    if(sum(reaction_species[,j][which(rownames(reaction_species) %in% temp_total_reactions[,1])]) == nrow(temp_total_reactions)){ # check if all reactions have a "1" in reaction species based on the set
      temp_species <- rbind(temp_species,colnames(reaction_species)[j]) # add species colname to temp_species table
    }else{
      temp_species <- temp_species
    }
    #temp_species <- na.omit(temp_species) # remove NA rows (which were added from the initial temp_species matrix)
    complete_pathways[[i]] <- temp_species
  }
  }else{
   ### constraint about transporter
    
    temp_species <- matrix(ncol=1) #temp file to store species that have complete pathway
    temp_total_reactions <- as.matrix(extract_total_reactions[[i]]) # number of reactions for each pathway
    transporters <- strsplit(extract_transporters[[i]],split=":") # split transporter cell by ":"
    
    for(tp in 1:length(transporters)){
      temp_total_reactions <- rbind(temp_total_reactions,transporters[[tp]][1]) # add reactions related to transporters to eaction list
    }
    
    
    for(j in 1:ncol(reaction_species)){ # iterate through each species in reaction_species table
      if(sum(reaction_species[,j][which(rownames(reaction_species) %in% temp_total_reactions[,1])]) == nrow(temp_total_reactions)){ # check if all reactions have a "1" in reaction species based on the set
        temp_species <- rbind(temp_species,colnames(reaction_species)[j]) # add species colname to temp_species table
      }else{
        temp_species <- temp_species
      }
      #temp_species <- na.omit(temp_species) # remove NA rows (which were added from the initial temp_species matrix)
      complete_pathways[[i]] <- temp_species
    }
  }
  
}

names(complete_pathways) <- list_names # Assign pathway names to each list

for(list in 1:length(complete_pathways)){
  ind <- apply(complete_pathways[[list]],1,function(x) all(is.na(x)))
  complete_pathways[[list]] <- complete_pathways[[list]][!ind,]
}

sink("./complete_pathway_species.txt") # save to file the species that have the complete ability to perform each pathway in user_defined data
print(complete_pathways)
sink() # close file


##### Script to obtain species capable of performing subsets of reactions ##########

complete_subsets <- list() # empty list to store results



for(list in 1:length(complete_pathways)){ # iterate thorugh each pathway in complete_pathways
  '
  if(length(complete_pathways[[list]])!=0){ # if there are species with full potential for a pathway
  temp_complete <- reaction_species[,-which(colnames(reaction_species) %in% complete_pathways[[list]])] # remove species with complete potential from subsequent steps
  }else{
    temp_complete <- reaction_species
  }
  '
  temp_complete <- reaction_species
  
  subsets_reactions <- strsplit(extract_subset_single_org[[list]],",") # split strings to get reactions for each subset
  
  if(length(subsets_reactions) != 0){ # check if there are any constraints for single orgs...
    
    temp_subset_list <- list()
    subset_names <- c()
    
    for(subset in 1:length(subsets_reactions)){ #iterate through each subset of reactions
    
      
    temp_species <- matrix(ncol=1) #temp file to store species that have complete subset of reactions
    
    temp_subset_reactions <- as.matrix(subsets_reactions[[subset]]) # temp file to store reactions for complete pathway
    subset_names <- append(subset_names,as.character(paste(temp_subset_reactions,collapse = ","))) # temp file to store reaction subset
    
    # Check if reactions in subsets require transporters
    
    temp_transp <- strsplit(unlist(extract_transporters[[list]]),split = ":|,")
    
    
    for(k in 1:length(temp_transp)){
      if(any(temp_subset_reactions %in% temp_transp[[k]])){
        temp_subset_reactions <- rbind(temp_subset_reactions, temp_transp[[k]][1])
      }else{
        temp_subset_reactions
      }
    }
  
       
        
      
    
    for(j in 1:ncol(temp_complete)){ # iterate through each species in temp_complete table
      if(sum(temp_complete[,j][which(rownames(temp_complete) %in% temp_subset_reactions[,1])]) == nrow(temp_subset_reactions)){ # check if all subset reactions have a "1" in reaction species based on the set
        temp_species <- rbind(temp_species,colnames(temp_complete)[j]) # add species colname to temp_complete table
      }else{
        temp_species <- temp_species
      }
     
      
    #temp_species <- temp_species[!is.na(temp_species),] # remove NA rows (which were added from the initial temp_complete matrix)
    temp_subset_list[[subset]] <- temp_species
    }
  }
    names(temp_subset_list) <- subset_names #Assign names of lists from reaction sets
  complete_subsets[[list]] <- temp_subset_list # add to list results
  complete_subsets[[list]] <- lapply(complete_subsets[[list]], function(x) x[!is.na(x)])
  
  }else{
  
  complete_subsets[[list]] <- NA
  }
}
names(complete_subsets) <- names(complete_pathways) # Assign the names of pathways to lists





'
Output of complete_subsets
==========================

$`Anaerobic benzoate-acetylCoA`
$`Anaerobic benzoate-acetylCoA`[[1]]
[1] "adv.faa" "azd.faa"

$`Anaerobic benzoate-acetylCoA`[[2]]
logical(0)


$`Aerobic benzoate-acetylCoA`
[1] NA
'

#######################################################################################
# Get potential interactions of species using subsets of reactions defined by the user including transporter constraints

single_org_interactions <- list() # Create empty list to store results

for(i in 1:length(complete_subsets)){ # for each pathway
  species_interactions <- list() # create temp file to store species with ability to perform reaction subsets
  
  
  check_lengths <- lengths(complete_subsets[[i]]) # Check if there were any subsets defined by the user that did not return any species
  
  if(all(check_lengths!=0) ==FALSE){
    single_org_interactions[[i]] <- paste("No microbial interactions")
  }else{
  for(j in 1:length(complete_subsets[[i]])){
    x <- unlist(complete_subsets[[i]][j]) # separate species in each reaction set
    species_interactions[[j]] <- x
    
  }
  single_org_interactions[[i]] <- expand.grid(species_interactions) # combine all species from all reaction sets
  colnames(single_org_interactions[[i]]) <- names(complete_subsets[[i]])
  }
  
}
names(single_org_interactions) <- names(complete_subsets) 

sink("./single_org_interactions.txt") # create file to store single_org_interactions
print(single_org_interactions)
sink() # close file


###################################################################################


# Transform table with output to json format files (path:reaction:enzyme)

#File required: user_input + gpr.xlsx
csv2 <- user_input[,c(1:2)]

temp_csv <- matrix(ncol = 1)
test_csv_path <- c()

for(l in 1:length(csv2)){
  z <- unlist(strsplit(toString(csv2[[2]][l]),split = ","))
  for(ll in 1:length(z)){
    temp_csv <- rbind(temp_csv,z[ll])
    test_csv_path <- append(test_csv_path,toString(csv2[[1]][l]),after = length(test_csv_path))
  }
  
}

ind <- apply(temp_csv, 1, function(x) is.na(x))
temp_csv <- as.matrix(temp_csv[ !ind, ])

test_csv <- cbind(test_csv_path,temp_csv,NA)

for(line in 1:nrow(test_csv)){
  test_csv[line,3] <- reac_ec_table[rownames(reac_ec_table) == test_csv[line,2],]
}

test_csv <- as.data.frame(test_csv)
colnames(test_csv) <- c("path","reaction","EC_number")

paths_final <- c()

for(path in 1:length(unique(test_csv$path))){
  path_ids <- unique(test_csv$path)
  x <- test_csv[test_csv$path %in% path_ids[path],]
  paths_file=c()
  reaction_ids <- unique(x$reaction)
  
  for(reaction in 1:length(unique(x$reaction))){
  
  y <- x[x$reaction %in% reaction_ids[reaction],]
  
  temp_paths <- paste("\"",x$reaction[reaction],"\": ","[","\"",y$EC_number,"\"", "]",sep="")
  
  paths_file <- append(paths_file,temp_paths,after = length(paths_file))
  }
  paths_file <- paste(paths_file,collapse = ",")
  total_paths <- paste("\"",x$path[path],"\": {",paths_file,"}",sep="")
  paths_final <- append(paths_final,total_paths,after = length(paths_final))
  
  paths_final <- paste(paths_final,collapse=",")
}
paths_final <- paste0("{",paths_final,"}")


sink("./paths.json")
cat(paths_final)
sink()

#########################################################################

# Transform table with output to json format files (enzyme:kos)

#input_file: ec_table


ec_file=c()
for(ec in 1:length(unique(rownames(ec_table)))){
  ec_ids <- unique(rownames(ec_table))
  x <- ec_table[rownames(ec_table) %in% ec_ids[ec],]
  
  ko_ids <- gsub(" ","",x)
  ko_ids <- strsplit(ko_ids,split=",")
  
  
  if(grepl("\\[", ko_ids[1]) == 0){ # Check if there is a "[" which indicates the "and" rule
    if(length(ko_ids[[1]])>1){
    ko_ids <- paste(ko_ids, sep=",")
    
    ko_ids <-gsub(pattern = "c\\(",replacement = "",x = ko_ids)
    ko_ids <-gsub(pattern = "\\)",replacement = "",x = ko_ids)
    temp_kos <- paste0("\"",ec_ids[ec], "\": ", "[", ko_ids, "]")
    }else{
      temp_kos <- paste0("\"",ec_ids[ec], "\": ", "[", "\"",ko_ids,"\"", "]")
    }
  }else{
    ko_ids <-gsub(pattern = "\\[",replacement = "",x = ko_ids)
    ko_ids <-gsub(pattern = "\\]",replacement = "",x = ko_ids)
    ko_ids <-gsub(pattern = "c\\(",replacement = "",x = ko_ids)
    ko_ids <-gsub(pattern = "\\)",replacement = "",x = ko_ids)
    
    temp_kos <- paste0("\"",ec_ids[ec], "\": ", "[[", ko_ids, "]]")
   
    
    
  }
  ec_file <- append(ec_file,temp_kos,after = length(ec_file))
}
ec_file <- paste0(ec_file,collapse = ",")
ec_file <- paste0("{",ec_file,"}")

sink("./GP_rules.json")
cat(ec_file)
sink()

#######################################################################


