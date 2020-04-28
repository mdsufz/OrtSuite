#!/bin/bash

# $1 - file Species_Annotation.csv (with path - should be found in the output folder 'Results of OrtAn
# $2 - file GP_rules.json (generated from the Rscript)
# $3 - file path.json (generated from the Rscript)
# $4 - file species_exclude.json (generated from the Rscript)

python ./microbial_interactions/main.py $1 $2 $3 $4

