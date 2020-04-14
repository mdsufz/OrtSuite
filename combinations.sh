#!/bin/bash

# $1 - file where python script is located
# $2 - file Species_Annotation.csv (with path - should be found in the output folder 'Results of OrtAn
# $3 - file GP_rules.json (generated from the Rscript)
# $4 - file path.json (generated from the Rscript)
# $5 - file species_exclude.json (generated from the Rscript)

python $1 $2 $3 $4 $5

