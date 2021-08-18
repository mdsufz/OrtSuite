#!/bin/bash
 
# $1 - project full path (It needs to be inside the OrtSuite main folder. Ex: "~/OrtSuite/PROJECT_NAME/")
# $2 - the flag of the list (-m, -r, -e, or -k) 
# $3 - full path of the identifiers list
# $4 - full path of Ortsuite instalation
 
mkdir "$1"/database
 
printf "%s\n" " " "seq_download_started" " "  
 
download_kos -o "$1"/database "$2" "$3"
 
printf "%s\n" " " "seq_download_finished" " "
 
# Copy "get_gpr.sh" from the OrtSuite folder to the database folder.
#cp "$1"/../get_gpr.sh "$1"/database

# if it is not in the main folder, run the command below:
cp "$4"/examples/test_database/get_gpr.sh "$1"/database
cd "$1"/database
printf "%s\n" "Starting GPR collection" " "

sh get_gpr.sh KO_grp.txt "$4"/GPRKO.jar final_gpr.xlsx
 
printf "%s\n" " " "Finished GPR collection - Do manual check of file" " "
 
# Done

