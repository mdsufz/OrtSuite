#!/bin/bash

find $1 -name '*.fa' > $2 #Find and save all KO ids from directory

sed -E 's/(.fa)+$//;s/.\///' $2 > $2 # remove ./ from the beginning of each line and .fa and store in a new file

java -jar $3 $2 $4
