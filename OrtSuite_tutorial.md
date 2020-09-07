Guide to predict putative microbial interactions with OrtSuite
====

Briefly, OrtSuite consists of **1:** Generating the Ortholog Reaction-Association database (*ORAdb*), **2:** Clustering of Orthologs, **3:** Functional annotation of clusters of orthologs, and **4:** Identification of putative microbial interactions.


1: Generating the Ortholog Reaction-Association database
====

Download all sequences associated with the list of enzyme commission numbers:

First make a directory: 

>mkdir examples/test_database # Make sure to be in the OrtSuite folder!

Download all sequences associated with provided list of reactions and generate ORAdb:
>download_kos -o ~/examples/test_database -r /OrtSuite/examples/complete_bta_reactions.txt

![download_kos](https://github.com/mdsufz/OrtSuite/blob/master/download_kos.png)

Extract GPR rules from KEGG:

>sh ./examples/test_database/get_gpr.sh ./examples/test_database/KO_gpr.txt ~/OrtSuite/keggOrthologues.jar ./examples/test_database/final_gpr.xlsx

The spreadsheet file should look like the following:

![download_kos](https://github.com/mdsufz/OrtSuite/blob/master/GPR_file.png)


2: Clustering of Orthologs
====

Perform clustering of genome sequences of interest with OrthoFinder :

>orthofinder -f ~/examples/orfs/ -o ~/examples/clusters/ -og

![orthofinder_results](https://github.com/mdsufz/OrtSuite/blob/master/orthofinder_result_folder.png)


3: Functional annotation of clusters of orthologs
====


Define the variables for input:

>work_dir="~/examples/OrtAn_results/"

>database="~/examples/test_database/"

>orthof="~/examples/clusters/Results_Apr24/"

**Note:** In this example a directory *clusters* was created to where the results from Orthofinder were copied. This is optional. 

>new_db="~/test_OrtSuite/new_db/" (optional)

Create the necessary directories:
>mkdir work_dir

>mkdir new_db (optional)

Create project:
>create_project -out $work_dir -db $database
 
Perform relaxed search:
>relaxed_search -wd $work_dir -of $orthof -t 2 -ident 40

Perform restrictive search:
>restrictive_search -wd $work_dir -t 2 -ident 70

Assign function to sequences in clusters of orthologs:
>annotation -wd $work_dir

![ortAn_results](https://github.com/mdsufz/OrtSuite/blob/master/ortAn_results_folder.png)

4: Identification of putative microbial interactions
====

Extract all microbial interactions with complete functional potential and using th constraints defined in *test_user_input.csv*:
```bash

Rscript gpr_manipulation.R -n ~/examples/test_database/final_gpr.xlsx -s ~/examples/OrtAn_new/Results/Species_Annotation.csv -u ~/examples/OrtAn_Results/Results/test_user_input.csv -o ~/examples/

```

