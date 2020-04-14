#!/bin/bash

# $1 - file gpr.xlsx (with path)
# $2 - file Species_Annotation.csv (with path - should be found in the output folder 'Results of OrtAn
# $3 - file where user input for constraints is  located (e.g. user_input.csv)

Rscript ./gpr_manipulation.R $1 $2 $3


