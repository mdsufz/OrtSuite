#!/bin/bash

# $1 - flag (-p , -m , -n)
# $2 - file gpr.xlsx or pathway or module (with path)
# $3 - file Species_Annotation.csv (with path - should be found in the output folder 'Results of OrtAn
# $4 - file where user input for constraints is  located (e.g. user_input.csv)

Rscript gpr_manipulation.R $1 $2 -s $3 -u $4


