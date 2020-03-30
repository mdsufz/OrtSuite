Guide to predict putative microbial interactions with OrtSuite
====

Briefly, OrtSuite consists of **1:** Generating the Ortholog Reaction-Association database (*ORAdb*), **2:** Clustering of Orthologs, **3:** Functional annotation of clusters of orthologs, and **4:** Identification of putative microbial interactions.


1: Generating the Ortholog Reaction-Association database
====

Download all sequences associated with the list of enzyme commission numbers:

First make a directory: 

>mkdir test_Ortsuite/test_database

Download all sequences associated with provided list of reactions and generate ORAdb:
>download_kos -o ~/test_OrtSuite/test_database -r /OrtSuite/examples/rx.txt

![download_kos](https://github.com/msdsufz/OrtSuite/blob/master/download_kos.png)

Extract GPR rules from KEGG:

>sh ./examples/test_database/get_gpr.sh . kos_4_gpr.txt ~/OrtSuite/keggOrthologues.jar gpr.xlsx

![download_kos](https://github.com/msdsufz/OrtSuite/blob/master/GPR_file.png)


2: Clustering of Orthologs
====

Perform clustering of genome sequences of interest with OrthoFinder :

>cd OrthoFinder

>orthofinder -f ~/Documents/Test_genomes/ -o ~/test_OrtSuite/orthofinder_res/ -og

![orthofinder_results](https://github.com/msdsufz/OrtSuite/blob/master/orthofinder_result_folder.png)

3: Functional annotation of clusters of orthologs
====


Define the variables for input:

>work_dir="~/test_OrtSuite/ortan_results/"

>database="~/test_OrtSuite/test_database/"

>orthof="~/test_OrtSuite/orthofinder_res/Results_Mar17_6/"

**Note:** In this example a directory *orthofinder_res* was created to where the results from Orthofinder were copied. This is optional. 

>new_db="~/test_OrtSuite/new_db/"

Create the necessary directories:
>mkdir work_dir
>mkdir new_db

Create project:
>create_project -out $work_dir -db $database
 
Perform relaxed search:
>relaxed_search -wd $work_dir -of $orthof -t 2 -ident 50

Perform restrictive search:
>restrictive_search -wd $work_dir -t 2

Assign function to sequences in clusters of orthologs:
>annotation -wd $work_dir -ident 95 -ppos 99 -qc 90 -sc 90

![ortAn_results](https://github.com/msdsufz/OrtSuite/blob/master/ortAn_results_folder.png)

4: Identification of putative microbial interactions
====

Extract all microbial interactions with complete potential:
```bash

Still missing.

```
Add constraints to reduce search space:

```
Still missing
```
