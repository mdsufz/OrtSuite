#!/bin/bash

# $1 - file gpr_manipulation.R
# $2 - flag (-p , -m , -n)
# $3 - file gpr.xlsx or pathway or module (with path)
# $4 - file Species_Annotation.csv (with path - should be found in the output folder 'Results of OrtAn
# $5 - file where user input for constraints is  located (e.g. user_input.csv)

Rscript $1 $2 $3 -s $4 -u $5


