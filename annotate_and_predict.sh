#!/bin/bash
 
# $1 - project folder's full path (ex: ~/OrtSuite/examples/ )
# $2 - orthof results folder's full path
# $3 - user_input.csv file's full path
# $4 - user-defined reaction pairs full path (Optional)
 
mkdir "$1"/work_dir
printf "%s\n"
 
work_dir="$1"/work_dir # location where you want to store the results
database="$1"/database # location of the ORAdb (where the FASTA files are located)
orthof="$2" # location of the folders with the results from OrthoFinder
user_input="$3"
 
create_project -out $work_dir -db $database
printf "%s\n"
 
relaxed_search -wd $work_dir -of $orthof -t 2
printf "%s\n"
 
restrictive_search -wd $work_dir -t 2
printf "%s\n"
 
annotation -wd $work_dir 
printf "%s\n"
 
mkdir $work_dir/interactions
printf "%s\n"
 
Rscript "$1"/../gpr_manipulation.R -n $database/final_gpr.xlsx -s $work_dir/Results/Species_Annotation.csv -u $user_input -o $work_dir/interactions

Rscript --vanilla "$1"/../network_visNetwork.R $work_dir/Results/Reactions_mapped_to_species.txt $4

printf "%s\n" "Done" " "
 
#Done
