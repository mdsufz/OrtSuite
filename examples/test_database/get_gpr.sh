#!/bin/bash

ls -d *.fa > temp.txt #Find and save all KO ids from directory

sed -E 's/(.fa)+$//;s/.\///' temp.txt > $1 # remove ./ from the beginning of each line and .fa and store in a new file

java -jar $2 $1 $3
