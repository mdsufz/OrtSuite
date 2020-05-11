#!/usr/bin/env Rscript


library("readxl")
library(stringr)
library(dplyr)
library(gsubfn)
library(mgsub)
library(stringr)
library(reshape)
library(reticulate)
library(optparse)

#output_ortan <- read.csv("/home/leonorfe/Documents/Test_genomes/total_bta/OrtAn_results/Results/Species_Annotation.csv",header=T,sep=";",row.names = 1) # Read in output file from OrtAn
#gpr_file <- read_xlsx("~/Documents/Test_genomes/total_bta/ORAdb/final_gpr.xlsx") # read in gpr file generated from .sh script
#user_input <- read.csv("~/Documents/Test_genomes/user_input.csv",header = T) # read in user input file for constraints 

###### ARGUMENT PARSING ######

option_list = list(
  make_option(c("-p", "--pathway_list"),
              type="character",
              default=NULL, 
              help="file with pathway identifier [default=%default]",
              metavar="character"),
  make_option(c("-m", "--module_list"),
              type="character",
              default=NULL,
              help="file with module identifiers [default=%default]",
              metavar="character"),
  make_option(c("-n", "--gpr_file"),
              type="character",
              default=NULL,
              help="all reactions in gpr_file - DEFAULT: final_gpr.xlsx [default=%default]",
              metavar="character"),
  make_option(c("-s", "--output_ortan_file"),
              type="character",
              default=NULL,
              help="Species_Annotation.csv file [default=%default]",
              metavar="character"),
  make_option(c("-u", "--user_input"),
              type="character",
              default=NULL,
              help="file with user-defined constraints [default=%default]",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default=NULL,
              help="output folder of resulting files [default=%default]",
              metavar="character"))



opt_parser = OptionParser(
  usage = "usage: %prog [options]",
  option_list=option_list,
  prog = "gpr_manipulation.R"
);

opt = parse_args(opt_parser);


if (is.null(opt$gpr_file)){
  print_help(opt_parser)
  stop("Please define gpr file", call.=FALSE)
}



# Define variables from parsed arguments
p_file <- opt$pathway_list # pathway file path
m_file <- opt$module_list # module file path
final_gpr_file <- opt$gpr_file # gpr_file path
output_ortan_file <- opt$output_ortan_file # Species_Annotation_file
user_input_file <- opt$user_input # user defined constraints file
output_folder <- opt$output #output folder for json files and combinations


output_ortan <- read.table(file=as.character(output_ortan_file),header=T,sep=";",row.names = 1) 


gpr_file <- read_xlsx(final_gpr_file) 


user_input <- read.csv(file=as.character(user_input_file),header = T) 



#If user provides pathway list

pathway_func <- function(){
#p_file <- "/home/leonorfe/test_pathway_optparse.txt"

pathway_file <- read.delim(file = p_file, header = F, sep = "\t", stringsAsFactors = F,colClasses = "character")

colnames(pathway_file) <- c("pathway_id")

p <- matrix(ncol=11) # Number of columns of gpr_file
colnames(p) <- colnames(gpr_file)

for(i in 1:nrow(gpr_file)){
  path_ids <- unlist(strsplit(gpr_file[i,]$`Pathway IDs`,split = "|",fixed = T))
  if(any(unlist(pathway_file) %in% path_ids)){
    p <- rbind(p,gpr_file[i,])  
  }else{
    p <- p
  }
}
p <- p[-1,]


### Get file mapping Reaction ids to KOs

temp_reaction_table <- unique(p$Reaction) # Get unique reaction ids

final_reaction_table <- matrix(NA,nrow=length(temp_reaction_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_reaction_table) <- temp_reaction_table
colnames(final_reaction_table) <- "GPR"

# paste all possible gpr rules for each unique reaction id

for (j in 1:length(temp_reaction_table)){
  x <- p[p$Reaction %in% temp_reaction_table[j],]
  
  final_reaction_table[j,1] <- unique(x$`Parsed Definitions Merged by Module and Reaction`)
  
}

final_reaction_table <<- final_reaction_table

### Get file mapping EC numbers to KOs


temp_ec_table <- unique(p$Protein) # Get unique EC numbers

final_ec_table <- matrix(NA,nrow=length(temp_ec_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_ec_table) <- temp_ec_table
colnames(final_ec_table) <- "kos"

# paste all possible gpr rules for each unique EC number

for (ec in 1:length(temp_ec_table)){
  x <- p[p$Protein %in% temp_ec_table[ec],]
  
  final_ec_table[ec,] <-  unique(x$`Parsed Definitions Merged by Module and Reaction`)
}

final_ec_table <<- final_ec_table

### Get file mapping Reaction ids to EC numbers


temp_reac_ec_table <- unique(p$Reaction) # Get unique reaction ids

final_reac_ec_table <- matrix(NA,nrow=length(temp_reac_ec_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_reac_ec_table) <- temp_reac_ec_table
colnames(final_reac_ec_table) <- "EC_number"

# paste all possible EC numbers for each unique reaction id

for (reaction in 1:length(temp_reac_ec_table)){
  x <- p[p$Reaction %in% temp_reac_ec_table[reaction],]
  if(nrow(x)==1){
    final_reac_ec_table[reaction,] <- x$Protein
  }else{
    temp_x <- paste0(shQuote(x$Protein),collapse=",") 
    temp_x <- sub('^.(.*).$','\\1',temp_x)
    final_reac_ec_table[reaction,] <- temp_x
  }
}

final_reac_ec_table <<- final_reac_ec_table

#Output file ---> EC numbers associated with each reaction
'
      EC_number  
R05575 "1.1.1.35","1.1.1.211" 
R03026 "4.2.1.17" 
R01778 "1.1.1.211"
R05593 "3.7.1.21" 
R02601 "4.2.1.80" 
R02451 "1.3.7.8"  
R00750 "4.1.3.39" 
R00813 "1.3.1.25" 
R05597 "4.2.1.100"
R02604 "3.7.1.9"  
R00238 "2.3.1.16" 
R00816 "1.13.11.2"
R05581 "1.1.1.368" 
'
}
################################################################

#If user provides module list

module_func <- function(){

#m_file <- "/home/leonorfe/test_module_optparse.txt"


module_file <- read.delim(file = m_file, header = F, sep = "\t", stringsAsFactors = F,colClasses = "character")

colnames(module_file) <- c("module_id")

m <- matrix(ncol=11)
colnames(m) <- colnames(gpr_file)

for(i in 1:nrow(gpr_file)){
  module_ids <- unlist(strsplit(gpr_file[i,]$Module, split = "|",fixed = T))
  if(any(unlist(module_file) %in% module_ids)){
    m <- rbind(m,gpr_file[i,])  
  }else{
    m <- m
  }
}

m <- m[-1,]


### Get file mapping Reaction ids to KOs

temp_reaction_table <- unique(m$Reaction) # Get unique reaction ids

final_reaction_table <- matrix(NA,nrow=length(temp_reaction_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_reaction_table) <- temp_reaction_table
colnames(final_reaction_table) <- "GPR"

# paste all possible gpr rules for each unique reaction id

for (j in 1:length(temp_reaction_table)){
  x <- m[m$Reaction %in% temp_reaction_table[j],]
  
  final_reaction_table[j,1] <- unique(x$`Parsed Definitions Merged by Module and Reaction`)
  
}

final_reaction_table <<- final_reaction_table

### Get file mapping EC numbers to KOs


temp_ec_table <- unique(m$Protein) # Get unique EC numbers

final_ec_table <- matrix(NA,nrow=length(temp_ec_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_ec_table) <- temp_ec_table
colnames(final_ec_table) <- "kos"

# paste all possible gpr rules for each unique EC number

for (ec in 1:length(temp_ec_table)){
  x <- m[m$Protein %in% temp_ec_table[ec],]
  
  final_ec_table[ec,] <-  unique(x$`Parsed Definitions Merged by Module and Reaction`)
}

final_ec_table <<- final_ec_table

### Get file mapping Reaction ids to EC numbers


temp_reac_ec_table <- unique(m$Reaction) # Get unique reaction ids

final_reac_ec_table <- matrix(NA,nrow=length(temp_reac_ec_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_reac_ec_table) <- temp_reac_ec_table
colnames(final_reac_ec_table) <- "EC_number"

# paste all possible EC numbers for each unique reaction id

for (reaction in 1:length(temp_reac_ec_table)){
  x <- m[m$Reaction %in% temp_reac_ec_table[reaction],]
  if(nrow(x)==1){
    final_reac_ec_table[reaction,] <- x$Protein
  }else{
    temp_x <- paste0(shQuote(x$Protein),collapse=",") 
    temp_x <- sub('^.(.*).$','\\1',temp_x)
    final_reac_ec_table[reaction,] <- temp_x
  }
}
final_reac_ec_table <<- final_reac_ec_table
}

############################################################

# if no list is provided use all GPR rules by merging them

no_func <- function(){

temp_reaction_table <- unique(gpr_file$Reaction) # Get unique reaction ids

final_reaction_table <- matrix(NA,nrow=length(temp_reaction_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_reaction_table) <- temp_reaction_table
colnames(final_reaction_table) <- "GPR"
  
# paste all possible gpr rules for each unique reaction id

for (j in 1:length(temp_reaction_table)){
  x <- gpr_file[gpr_file$Reaction %in% temp_reaction_table[j],]
  
  final_reaction_table[j,1] <- unique(x$`Parsed Definitions Merged by Module and Reaction`)

}

final_reaction_table <<- final_reaction_table


# Output file ---> GPR rules for each reaction
'
       Definition            
R00238 "[[K00632, K07508, K07509, K07513],[K000626]]"                                                                                         
R00750 "[[K01666],[K18365, K01666]]"                                                                                                          
R00813 "[[K05783]]"                                                                                                                           
R00816 "[[K00446, K07104]]"                                                                                                                   
R01778 "[[K07511, K07515],[K07514, K07515, K07511],[K00022, K07516, K01825, K01782, K07514, K07515, K10527],[K07511, K07514, K07515, K14729]]"
R01976 "[[K00074]]"                                                                                                                           
R02451 "[[[K04112,K04113,K04114,K04115]],[[K19515,K19516]]]"                                                                                  
R02487 "[[K00252]]"                                                                                                                           
R02601 "[[K02554],[K18364, K02554]]"                                                                                                          
R02604 "[[K01617, K10216]]"                                                                                                                   
R03026 "[[K01692, K07511, K13767, K01825, K01782, K07514, K07515, K10527],[K01692]]"                                                          
R05575 "[[K07547,K07548]]"                                                                                                                    
R05581 "[[K07538]]"                                                                                                                           
R05593 "[[K07539]]"                                                                                                                           
R05597 "[[K07537]]"

'


############ Obtain kos associated with each EC number ################


temp_ec_table <- unique(gpr_file$Protein) # Get unique EC numbers

final_ec_table <- matrix(NA,nrow=length(temp_ec_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_ec_table) <- temp_ec_table
colnames(final_ec_table) <- "kos"

# paste all possible gpr rules for each unique EC number

for (ec in 1:length(temp_ec_table)){
  x <- gpr_file[gpr_file$Protein %in% temp_ec_table[ec],]
  
  final_ec_table[ec,] <-  unique(x$`Parsed Definitions Merged by Module and Reaction`)
}


final_ec_table <<- final_ec_table

#Output file  ---> GP rules for each enzyme (will be used when transforming to .json)
'

            kos                   
2.3.1.9   "[[K00632, K07508, K07509, K07513],[K000626]]"                                                                                         
2.3.1.16  "[[K00632, K07508, K07509, K07513],[K000626]]"                                                                                         
4.1.3.39  "[[K01666],[K18365, K01666]]"                                                                                                          
1.3.1.25  "[[K05783]]"                                                                                                                           
1.13.11.2 "[[K00446, K07104]]"                                                                                                                   
1.1.1.211 "[[K07511, K07515],[K07514, K07515, K07511],[K00022, K07516, K01825, K01782, K07514, K07515, K10527],[K07511, K07514, K07515, K14729]]"
1.1.1.157 "[[K00074]]"                                                                                                                           
1.3.7.8   "[[[K04112,K04113,K04114,K04115]],[[K19515,K19516]]]"                                                                                  
1.3.8.6   "[[K00252]]"                                                                                                                           
4.2.1.80  "[[K02554],[K18364, K02554]]"                                                                                                          
3.7.1.9   "[[K01617, K10216]]"                                                                                                                   
4.2.1.17  "[[K01692, K07511, K13767, K01825, K01782, K07514, K07515, K10527],[K01692]]"                                                          
1.1.1.35  "[[K07547,K07548]]"                                                                                                                    
1.1.1.368 "[[K07538]]"                                                                                                                           
3.7.1.21  "[[K07539]]"                                                                                                                           
4.2.1.100 "[[K07537]]"

'
############# Obtain reaction associated with EC numbers ##############



temp_reac_ec_table <- unique(gpr_file$Reaction) # Get unique reaction ids

final_reac_ec_table <- matrix(NA,nrow=length(temp_reac_ec_table)) # Create a matrix to store GPR rules for unique reaction ids
row.names(final_reac_ec_table) <- temp_reac_ec_table
colnames(final_reac_ec_table) <- "EC_number"

# paste all possible EC numbers for each unique reaction id

for (reaction in 1:length(temp_reac_ec_table)){
  x <- gpr_file[gpr_file$Reaction %in% temp_reac_ec_table[reaction],]
  
    if(length(unique(x$Protein))==1){
      final_reac_ec_table[reaction,] <- unique(x$Protein)
    }else{
      temp_x <- paste0(shQuote(unique(x$Protein)),collapse=",") 
      temp_x <- sub('^.(.*).$','\\1',temp_x)
      final_reac_ec_table[reaction,] <- temp_x
    }
  }


final_reac_ec_table <<- final_reac_ec_table

#Output file ---> EC numbers associated with each reaction
'
      EC_number  
R00238 "2.3.1.9;2.3.1.16"
R00750 "4.1.3.39"        
R00813 "1.3.1.25"        
R00816 "1.13.11.2"       
R01778 "1.1.1.211"       
R01976 "1.1.1.157"       
R02451 "1.3.7.8"         
R02487 "1.3.8.6"         
R02601 "4.2.1.80"        
R02604 "3.7.1.9"         
R03026 "4.2.1.17"        
R05575 "1.1.1.35"        
R05581 "1.1.1.368"       
R05593 "3.7.1.21"        
R05597 "4.2.1.100" 
'
}


####################################
#    check which option has flag
####################################


if(is.null(opt$p) & is.null(opt$m)){
  print("using complete set of reactions")
  no_func()
}else if(is.null(opt$p)){
  print("using module")
  module_func()
}else if(is.null(opt$m)){
  print("using pathway")
  pathway_func()
  }



###################################################################################
##         Create json files for Marta's tool                                   ###
###################################################################################

# Transform table with output to json format files (path:reaction:enzyme)

csv2 <- user_input[,c(1:2)]

temp_csv <- matrix(ncol = 1)
test_csv_path <- c()

for(l in 1:nrow(csv2)){
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
  
  if(length(final_reac_ec_table[rownames(final_reac_ec_table) == test_csv[line,2]]) == 0){
  test_csv[line,3] <- 0
  }else{
  test_csv[line,3] <- final_reac_ec_table[rownames(final_reac_ec_table) == test_csv[line,2],]
  }
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
paths_final <- gsub(pattern = "\'",replacement = "\"",paths_final)

sink(paste0(output_folder,"/paths.json",collapse = ""))
cat(paths_final)
sink()

#########################################################################

# Transform table with output to json format files (enzyme:kos)


ec_file=c()
for(ec in 1:length(unique(rownames(final_ec_table)))){
  ec_ids <- unique(rownames(final_ec_table))
  x <- final_ec_table[rownames(final_ec_table) %in% ec_ids[ec],]
  
  ko_ids <- gsub(" ","",x) 
  ko_ids <- sub('^.(.*).$','\\1',ko_ids) 
  
  ko_ids <- strsplit(ko_ids,split="\\],") 
  temp_ec_file <- c()
  
  if(length(ko_ids[[1]])>1){
   
     temp_ko_and <- c()
     
     
    for(i in 1:length(ko_ids[[1]])){
      if(grepl("\\[", ko_ids[[1]][i]) == 0){
        
        temp_ko_ids <- strsplit(ko_ids[[1]],split=",")
        temp_ko_ids <-gsub("\\[","",temp_ko_ids)
        temp_ko_ids <-gsub("\\]","",temp_ko_ids)
        
        temp_ko_ids <-gsub(pattern = "c\\(",replacement = "",x = temp_ko_ids)
        temp_ko_ids <-gsub(pattern = "\\)",replacement = "",x = temp_ko_ids)
        
        temp_ko_ids <- paste0("[", temp_ko_ids, "]")
        
      } else{
        if(grepl("\\[\\[", ko_ids[[1]][i])==0){
        temp_ko_ids <- strsplit(ko_ids[[1]][i],split=",")
        temp_ko_ids <-gsub("\\[","",temp_ko_ids)
        temp_ko_ids <-gsub("\\]","",temp_ko_ids)
        
        temp_ko_ids <-gsub(pattern = "c\\(",replacement = "",x = temp_ko_ids)
        temp_ko_ids <-gsub(pattern = "\\)",replacement = "",x = temp_ko_ids)
        
        temp_ko_ids <- paste0("[", temp_ko_ids, "]")
        }else{
          temp_ko_ids <- strsplit(ko_ids[[1]][i],split=",")
          temp_ko_ids <-gsub("\\[","",temp_ko_ids)
          temp_ko_ids <-gsub("\\]","",temp_ko_ids)
          
          temp_ko_ids <-gsub(pattern = "c\\(",replacement = "",x = temp_ko_ids)
          temp_ko_ids <-gsub(pattern = "\\)",replacement = "",x = temp_ko_ids)
          
          temp_ko_ids <- paste0("[[", temp_ko_ids, "]]")
        } 
      }
      temp_ko_and <- append(temp_ko_and,temp_ko_ids,after = length(temp_ko_and))
    }
     temp_ko_and <- paste0(temp_ko_and,collapse = ",")
     temp_ko_final <- paste0("\"",ec_ids[ec], "\": ","[",temp_ko_and,"]")
     
     ec_file <- append(ec_file,temp_ko_final,after = length(ec_file))
    
  }else{
    if(grepl("\\[", ko_ids[[1]]) == 0){
      temp_ko_ids <- strsplit(ko_ids[[1]],split=",")
      
      if(length(temp_ko_ids[[1]])>1){
        temp_ko_ids <-gsub("\\[","",temp_ko_ids)
        temp_ko_ids <-gsub("\\]","",temp_ko_ids)
        temp_ko_ids <-gsub(pattern = "c\\(",replacement = "",x = temp_ko_ids)
        temp_ko_ids <-gsub(pattern = "\\)",replacement = "",x = temp_ko_ids)
    
        temp_ko_ids <- paste0("[", temp_ko_ids, "]")
      }else{
        init_bracket <- "["
        end_bracket <- "]"
        temp_ko_ids <- paste(init_bracket,temp_ko_ids,end_bracket,sep = "\"")
      }
    }else{
      temp_ko_ids <- strsplit(ko_ids[[1]],split=",")
      temp_ko_ids <-gsub("\\[","",temp_ko_ids)
      temp_ko_ids <-gsub("\\]","",temp_ko_ids)
      
      temp_ko_ids <-gsub(pattern = "c\\(",replacement = "",x = temp_ko_ids)
      temp_ko_ids <-gsub(pattern = "\\)",replacement = "",x = temp_ko_ids)
      
      temp_ko_ids <- paste0("[[", temp_ko_ids, "]]")
  }
  temp_ko_final <- paste0("\"",ec_ids[ec], "\": ",temp_ko_ids)
  ec_file <- append(ec_file,temp_ko_final,after = length(ec_file))
}




  
  }
ec_file <- paste0(ec_file,collapse = ",")
ec_file <- paste0("{",ec_file,"}")

  
  
sink(paste0(output_folder,"/GP_rules.json",collapse = ""))
cat(ec_file)
sink()

#######################################################################


########## Create table mapping reactions to species based on gpr rules  ###############


reaction_species <- matrix(NA,nrow = nrow(final_reaction_table),ncol=ncol(output_ortan))

rownames(reaction_species) <- rownames(final_reaction_table)
colnames(reaction_species) <- colnames(output_ortan)

for(i in 1:nrow(reaction_species)){
  kos <- final_reaction_table[i,1] 
  
  kos <- sub('^.(.*).$','\\1',kos) 
  x <- strsplit(kos,",\\[") 
  
  for(j in 1:ncol(output_ortan)){ 
    species_results <- as.matrix(as.data.frame(output_ortan[,j],row.names = rownames(output_ortan))) 
    
      if(grepl("\\[",x)==T){ 
        temp_x <- strsplit(x[[1]],",")
        for(list in 1:length(temp_x)){
          if(grepl("\\[",temp_x[list])==F){
            z <- as.character(strsplit(temp_x[[list]]," "))
            if(sum(species_results[which(rownames(species_results) %in% z)]) >= 1){ 
              reaction_species[i,j] <- 1
            }else{
              reaction_species[i,j] <- 0
            }
          }else{
            z <- as.character(strsplit(temp_x[[list]]," "))
            z <- gsub("\\[","",z) 
            z <- gsub("\\]","",z) 
            if(sum(species_results[which(rownames(species_results) %in% z)]) == length(z)){ 
               reaction_species[i,j] <- 1 
               break
            }else{
               reaction_species[i,j] <- 0
            }
          }
        }
    }else{
      kos <- gsub(" ","",kos) 
      x <- strsplit(kos,",") 
      if(sum(species_results[which(rownames(species_results) %in% x[[1]])]) >= 1){ 
        reaction_species[i,j] <- 1
      }else{
        reaction_species[i,j] <- 0
      }
    }
  }
}


write.csv(reaction_species,file=paste0(output_folder,"/Reactions_mapped_to_species.csv",collapse = ""))



#### Table with user defined constraints ########
  

extract_total_reactions <- strsplit(as.character(user_input$Reactions),",") 
extract_subset_single_org <- str_match_all(user_input$Single_org, "(?<=\\().+?(?=\\))") 
extract_transporters <- str_match_all(user_input$Transporter, "(?<=\\().+?(?=\\))") 
names(extract_transporters) <- user_input$Pathway 

#### Script to obtain all species individually capable of performing complete pathways ##########

complete_pathways <- list() 


list_names <- as.character(user_input$Pathway) 

for(i in 1:length(extract_total_reactions)){ 
  if(length(extract_transporters[[i]])==0){
    
  temp_species <- matrix(ncol=1) 
  
  temp_total_reactions <- as.matrix(extract_total_reactions[[i]]) 
  
  for(j in 1:ncol(reaction_species)){
    if(sum(reaction_species[,j][which(rownames(reaction_species) %in% temp_total_reactions[,1])]) == nrow(temp_total_reactions)){ 
      temp_species <- rbind(temp_species,colnames(reaction_species)[j]) 
    }else{
      temp_species <- temp_species
    }
    
    complete_pathways[[i]] <- temp_species
  }
  }else{
  
    
    temp_species <- matrix(ncol=1) 
    temp_total_reactions <- as.matrix(extract_total_reactions[[i]]) 
    transporters <- strsplit(extract_transporters[[i]],split=":") 
    
    for(tp in 1:length(transporters)){
      temp_total_reactions <- rbind(temp_total_reactions,transporters[[tp]][1]) 
    }
    
    ind_na <- apply(temp_total_reactions,1,function(x) all(is.na(x)))
    temp_total_reactions <- as.matrix(temp_total_reactions[!ind_na,])
    
    for(j in 1:ncol(reaction_species)){ 
      if(sum(reaction_species[,j][which(rownames(reaction_species) %in% temp_total_reactions[,1])]) == nrow(temp_total_reactions)){ 
        temp_species <- rbind(temp_species,colnames(reaction_species)[j]) 
      }else{
        temp_species <- temp_species
      }
    
      complete_pathways[[i]] <- temp_species
    }
  }
  
}

names(complete_pathways) <- list_names 

for(list in 1:length(complete_pathways)){
  ind <- apply(complete_pathways[[list]],1,function(x) all(is.na(x)))
  complete_pathways[[list]] <- complete_pathways[[list]][!ind,]
}

sink(paste0(output_folder,"/complete_pathway_species.txt",collapse = "")) 
print(complete_pathways)
sink() 

#### Get all potential microbial interactions (excluding those with the ability to perform complete pathways) based solely on the reaction presence

for(list in 1:length(complete_pathways)){ 
  
  if(length(complete_pathways[[list]])!=0){ 
    complete_species <- reaction_species[,-which(colnames(reaction_species) %in% complete_pathways[[list]])] 
  }else{
    complete_species <- reaction_species
  }
  
  ## subset table to only include the reactions in pathway
  complete_species <- as.matrix(complete_species[rownames(complete_species) %in% extract_total_reactions[[list]],])
  
  # First, get all combinations of column names, with varying lengths (vary from 1 to 18 columns), and put into a list.
  lst <- do.call(c, lapply(seq_along(colnames(complete_species)), combn, x = colnames(complete_species), simplify = FALSE))
 
  # Second, would go through that list, and create matrices for those selected columns.
  lst2 <- lapply(lst, function(x) matrix(complete_species[,x], ncol = length(x)))
  
  # Third, would name the list based on the column names and remove the quotes.
  names(lst2) <- gsub('[\"]', '', lst)
  
  # Fourth, would check each matrix in the list where (a) any column in a row has a 1, and (b) all rows meet criteria (a).
  reaction_comb_species <- names(which(lapply(lst2, function(x) { all(apply(x, 1, function(y) any(y > 0))) } ) == TRUE))
  
  # Write combinations that meet the criteria to file.
  write.csv(reaction_comb_species,file=paste0(output_folder,"/",names(complete_pathways[list]),"_interactions.csv",collapse = ""))
   
}
  


##### Script to obtain species capable of performing subsets of reactions ##########

complete_subsets <- list() 



for(list in 1:length(complete_pathways)){ 
  
  if(length(complete_pathways[[list]])!=0){ 
  temp_complete <- reaction_species[,-which(colnames(reaction_species) %in% complete_pathways[[list]])] 
  }else{
    temp_complete <- reaction_species
  }
  
  
  
  subsets_reactions <- strsplit(extract_subset_single_org[[list]],",") 
  
  if(length(subsets_reactions) != 0){ 
    
    temp_subset_list <- list()
    subset_names <- c()
    
    for(subset in 1:length(subsets_reactions)){ 
    
      
    temp_species <- matrix(ncol=1) 
    
    temp_subset_reactions <- as.matrix(subsets_reactions[[subset]]) 
    subset_names <- append(subset_names,as.character(paste(temp_subset_reactions,collapse = ","))) 
    
    
    temp_transp <- strsplit(unlist(extract_transporters[[list]]),split = ":|,")
    
    
    for(k in 1:length(temp_transp)){
      if(any(temp_subset_reactions %in% temp_transp[[k]])){
        temp_subset_reactions <- rbind(temp_subset_reactions, temp_transp[[k]][1])
      }else{
        temp_subset_reactions
      }
    }
  
     
    for(j in 1:ncol(temp_complete)){ 
      if(sum(temp_complete[,j][which(rownames(temp_complete) %in% temp_subset_reactions[,1])]) == nrow(temp_subset_reactions)){ 
        temp_species <- rbind(temp_species,colnames(temp_complete)[j]) 
      }else{
        temp_species <- temp_species
      }
     
      
    temp_subset_list[[subset]] <- temp_species
    }
  }
    names(temp_subset_list) <- subset_names 
  complete_subsets[[list]] <- temp_subset_list 
  complete_subsets[[list]] <- lapply(complete_subsets[[list]], function(x) x[!is.na(x)])
  
  }else{
  
  complete_subsets[[list]] <- NA
  }
}
names(complete_subsets) <- names(complete_pathways) 



#######################################################################################


single_org_interactions <- list() 

for(i in 1:length(complete_subsets)){ 
  species_interactions <- list() 
  
  
  check_lengths <- lengths(complete_subsets[[i]]) 
  
  if(all(check_lengths!=0) ==FALSE){
    single_org_interactions[[i]] <- paste("No microbial interactions possible under defined contraints")
  }else{
  for(j in 1:length(complete_subsets[[i]])){
    x <- unlist(complete_subsets[[i]][j])
    species_interactions[[j]] <- x
    
  }
  single_org_interactions[[i]] <- expand.grid(species_interactions) 
  colnames(single_org_interactions[[i]]) <- names(complete_subsets[[i]])
  }
  
}
names(single_org_interactions) <- names(complete_subsets) 

sink(paste0(output_folder,"/single_org_subset_interactions.txt",collapse = "")) 
print(single_org_interactions)
sink()

####


