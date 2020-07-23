# OrtSuite - a flexible pipeline for annotation of ecosystem processes and prediction of putative microbial interactions

OrtSuite was developed with the goal to facilitate annotation of ecosystem processes and identify putative microbial interactions by automating the complete process from sequence retrieval, clustering of ortholog sequences, functional annotation, to putative microbial interactions prediction. 



![workflow](https://github.com/mdsufz/OrtSuite/blob/master/workflow_ortSuite.png)

**OrtSuite workflow** 

**(a)** Protein sequence from samples supplied by the user are clustered using OrthoFinder. 

**(b)** OrtScraper takes a text file containing a list of identifiers for each reaction in the pathway of interest supplied by the user to retrieve all protein sequences from KEGG.

**(c)** Sequences mapped to reactions are stored in ORAdb.

**(d)** Functional annotation consists of a two-stage process (relaxed and restrictive search). Relaxed search **(e)** performs sequence alignments between 10% of randomly selected sequences from each generated cluster. Clusters whose representative sequences share a minimum 50% identity to sequences in reaction set(s) in ORAdb transition to the restrictive search **(f)**. Here, all sequences from the cluster is aligned to all sequences in the corresponding reaction set(s) to which they had a hit. 

**(g and h)** The annotated sequences are used to identify putative microbial interactions based on their functional potential.

**(rounded blue rectangles)** Additional constraints can be added to reduce the search space of microbial interactions.


# Overview of OrtSuite

**Installation**

**ORAdb construction:** Generation of the user-defined Ortholog-Reaction Association (ORAdb) database

**OrthoFinder:** Clustering of orthologs  

**OrtScraper:** Bulk download of protein sequences for populating a user-defined database  

**OrtAn:** Functional annotation of clusters of orthologs  

**Interspecies interactions:** Prediction of interspecies interactions based on the functional potential of individual species  


# System Requirements

Resources for Ortsuite will vary depending on the amount of data being processed. In the example provided (consisting of 7 reactions with 16 associated KEGG Ortholog identifiers), we used an Intel Core i5-6200U 2.3GHz with 4 cores and 16 Gb of RAM. OrtSuite officially supports only Linux OS. 

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Installation
============


# OrtSuite 

OrtSuite is a python tool that performs functional annotation of clusters of orthologs and identifies putative microbial interactions. This tool automatically retrieves sequence data in bulk from KEGG (Kyoto Encyclopedia of Genes and Genomes) database to generate the Ortholog Reaction-Associated user-defined database (*ORAdb*).
Generation of clusters of orthologs is performed by OrthoFinder.

**Requirements:**  Python 3.6

**Dependencies:**  bs4, grequests, [OrthoFinder](https://github.com/davidemms/OrthoFinder), [DIAMOND](https://github.com/bbuchfink/diamond)


We suggest the use of the package manager [pip](https://pip.pypa.io/en/stable/) to install the virtual environment (*virtualenv*).


```bash
python3.6 -m pip install virtualenv
```

Create a virtual environment.

```bash
virtualenv venv_OrtSuite
```

Activate the virtual environment.

```bash
source venv_OrtSuite/bin/activate
```

Install dependencies

```bash
pip install grequests
pip install bs4
```

Next move to the folder where the file setup.py from the OrtSuite tool is located.

```bash
cd /path/to/OrtSuite
```

Run the command:

```bash
python3.6 setup.py install
```

Once the installation is finished the tool should be ready to use.

**Note:** Some users might encounter an error related to SSL certificates. If this is the case then type the following in your terminal:

> sudo update-ca-certificates --fresh

> export SSL_CERT_DIR=/etc/ssl/certs


OrthoFinder
=====

Install OrthoFinder (Available here: https://github.com/davidemms/OrthoFinder)
```bash
Download latest version: wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz 
Decompress file: tar xzf OrthoFinder.tar.gz (to your desired directory)
Testing: ./OrthoFinder/orthofinder -h (help menu)
```
Dependencies of OrthoFinder

MCL
====

The mcl clustering algorithm is available in the repositories and can be installed as any other package.
```bash
sudo apt-get install mcl
```



DIAMOND
====

Install DIAMOND (Available here: https://github.com/bbuchfink/diamond/releases)
```bash
Download latest version: wget https://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz
Decompress file: tar xzf diamond-linux64.tar.gz
Copy DIAMOND to your local directory: sudo cp diamond /usr/local/bin
```

BLAST+
====

Orthofinder allows the user to use BLAST+ instead of DIAMOND.
To install BLAST+ use:

```bash
sudo apt-get install ncbi-blast+
```

*Note* - No other dependency of OrthoFinder is required to run OrtSuite1.0 since we only use it to perform the clustering of orthologs. 


JAVA
====

To download the Gene-Protein-Reaction rules you need to have java installed. 

```bash
sudo apt install default-jre
```
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Using OrtSuite
=====

***For viewing an example of how to use OrtSuite please click [tutorial](OrtSuite_tutorial.md).***

Once installation of OrtSuite and all dependencies are completed the different commands can be called independently.



## Step 1) Generating the Ortholog Reaction-Assosiation database (ORAdb)


The generation of a user-defined ortholog reaction-association database starts with the download of all the sequences associated with the KO (KEGG Orthology) group to a FASTA file.

The input can be:

- a KEGG pathway map ID (e.g. map00362)

- [List of KO IDs](https://github.com/mdsufz/OrtSuite/blob/master/examples/kos.txt)

- [List of KEGG Reaction IDs](https://github.com/mdsufz/OrtSuite/blob/master/examples/rx.txt)

- [List of EC (Enzyme commission) numbers](https://github.com/mdsufz/OrtSuite/blob/master/examples/ecs.txt)

Note: The format of the input lists must be a txt file with only one ID per line. 

Run the command with the help option -h to see the usage and all the available options.

```bash
download_kos -h

usage: download_kos [-h] [--version] -o OUTPUT_DIR [-s SIZE]
                    (-m MAP | -r REACTIONS | -e EC_NUMBERS | -k KOS) [-p | -g]
                    [-v]

Create a database with all the KO (KEGG Orthology) sequences associated with
the give ID/s (KEGG pathway ID, KEGG Rections IDs, EC numers or KEGG orthology
IDs)

optional arguments:
  -h, --help     show this help message and exit
  --version      show program's version number and exit
  -o OUTPUT_DIR  Output directory to store the database.
  -s SIZE        Size - defines the number of requests that are made to the
                 KEGG database at the same time. DEFAULT: 5
  -m MAP         KEGG pathway ID
  -r REACTIONS   Path to txt file containing KEGG reaction IDs (one per each
                 line)
  -e EC_NUMBERS  Path to txt file containing EC numbers (one per each line)
  -k KOS         Path to txt file containing KEGG Orthology IDs (one per each
                 line)
  -p             Use this option if you want to download amino acid sequences.
                 (DEFAULT)
  -g             Use this option if you want to download nucleotide sequences.
  -v, --verbose  set loglevel to DEBUG
```

Examples of commands: 
====
When using a single pathway map from KEGG
```bash
download_kos -o output_folder -m map00362
```
When using reaction identifiers
```bash
download_kos -o output_folder -r reaction_list.txt 
```
When using EC numbers
```bash
download_kos -o output_folder -e ec_list.txt
```
When using KO identifiers
```bash
download_kos -o output_folder -k ko_list.txt
```
**Be aware that the output_folder must be previously created by the user!**

To test if the tool is working you can use the files contained in the [examples](https://github.com/mdsufz/OrtSuite/blob/master/examples) folder.

Note: Running this command may take some time and memory space.


ORAdb Output
======
A FASTA format file for each selected KO will be placed in the output folder defined by the user.
If you use -e or -r option, an additional file will be generated *associations.txt*, which indicates which kos where selected for download for each reaction/EC number.
In the same folder a file *info_db.csv* contains a table with information regarding the KO's that were selected for download along with their name and associated EC numbers and Reactions IDs. 
In the case of using KO identifiers the file associations.txt must be manually added (for an example please see [associations.txt](examples/associations.txt)).

Step 2) Downloading Gene-Protein-Reaction (GPR) rules
====

Once the FASTA files containing the sequences for the list of KOs is completed you can obtain the GPR rules.
First, copy the *get_gpr.sh* to the ORAdb output_folder you previously created and run the following command:

```bash
sh get_gpr.sh path/to/folder/with/KO_gpr.txt /path/to/path/to/keggOrthologues.jar /path/to/folder/for/final_gpr.xlsx
```

The get_pr.sh script reads all *.fa* files in the output_folder and calls an additional java script ([keggOrthologues.jar](keggOrthologues.jar)) which retrieves from KEGG Modules the respective KO-related GPR rules. *KO_gpr.txt* is a text file created by the script to store a list of all KOs in database and *final_gpr.xlsx* is a file created with the GPR rules.  

**Note 1:** Both *KO_gpr.txt* and *final_gpr.xlsx* are suggested filenames and can be altered by the user.  
**Note 2:** The *final_gpr.xlsx* file needs to have the *xlsx* extension.  
**Note 3:** Please be aware that the [cache folder](cache) and [cache.ccf](cache.ccf) file needs to stored in the same folder as [keggOrthologues.jar](keggOrthologues.jar). 

**Due to the limited information in KEGG database concerning GPR rules, manual inspection of the final_gpr.xlsx is advised! For a more comprehensive explanation please see the tutorial with an example.**  

Step 3) OrthoFinder
====

OrthoFinder takes as input a folder containing the FASTA sequences the user wants to cluster.

```bash
~/OrthoFinder/orthofinder -f ~/path/to/sequence/folder -o /path/to/output/folder -og
```

Note: If you wish to use BLAST+ instead of DIAMOND please use the following:

```bash
~/OrthoFinder/orthofinder -f ~/path/to/sequence/folder -o /path/to/output/folder -S blast -og
```
**Note:** OrthoFinder's output folder is generated automatically. The user can, however, define the parent directory where to store the output folder (e.g. /Documents/).

Step 4) OrtAn
====

OrtAn uses the output from [OrthoFinder](https://github.com/davidemms/OrthoFinder) to perform the annotation of the generated clusters of orthologs based on the user-defined ORAdb.

Overview of OrtAn:

>Create Project: OrtAn receives the information of the ORAdb database
>               and creates the working directory structure necessary for the subsequent tasks.

>Relaxed Search: This task takes the OrthoFinder information and runs 
>               a first relaxed search with DIAMOND to identify the associations 
>               between the returned clusters of orthologs and the different functions of the ORA database.

>Restrictive Search: This task performs a restrictive search only between 
>               the clusters of orthologs and groups of functions from the ORA 
>               database that were related during the relaxed_search.

>Annotation: This task consists in the annotation of the sequences present
>           in the clusters of orthologs after performing the restrictive_search.
>           Additionally the user has the possibility to create a new database with
>           the results or update the user-defined ORAdb (*create_db*).

>Identification of putative microbial interactions: This task allows the user to extract
>           the putative microbial interactions based on different sets of constraints 
>           (e.g. number of interacting species, ability to perform complete pathways, 
>           groups of interactions where a single species is responsible for a subset of reactions in the pathway).

### Inputs:

- OrthoFinder output directories    
- ORAdb  


Before running OrtAn you need to:  

- Run OrthoFinder with the input genomes  
- Prepare the database with the necessary format.  

### Required files and directories from OrthoFinder output:

- [Orthogroups.txt](examples/Orthogroups/Orthogroups.txt)  ```(located in examples/Orthogroups/Orthogroups.txt)```

- [SpeciesIDs.txt](examples/WorkingDirectory/SpeciesIDs.txt)  ```(located in examples/WorkingDirectory/SpeciesIDs.txt)```

- [SequenceIDs.txt](examples/WorkingDirectory/SequenceIDs.txt)  ```(located in examples/WorkingDirectory/SequenceIDs.txt)```

- [Orthogroup_Sequences](examples/Orthogroup_Sequences) folder  ```(located in examples/Orthogroup_Sequences/)```


## Create Project

The first task in OrtAn is to create the folder where all the results will be stored. During this task the location of the ORAdb is also defined. 

**Note:** you should indicate always the same output/working directory in all the steps.

Run ```create_project -h``` to see the usage of this command.

```bash
usage: create_project [-h] -out OUTPUT -db DATA [-l] [-v]

Creates a new project and returns the working directory to use in the other
steps of the pipline. (Necessary to access the temporary data that the tool
stores)

optional arguments:
  -h, --help            show this help message and exit
  -out OUTPUT, --output OUTPUT
                        Path to folder to store the output results
  -db DATA, --data DATA
                        Path to folder containing the fasta files with
                        sequences from the functions of interest
  -l, --logfile         To send log messages to a file int the output
                        directory
  -v, --verbose         set loglevel to DEBUG
```


## Relaxed Search

The first annotation of sequences in the clusters (*relaxed search*) only requires the user to define the working folder (also used in *create_project*) and the location of the folders genertaed by OrthoFinder.  

Ten (10) percent of random sequences from each cluster is used to retrieve the associations which will reduce the search space on the second annotation *restrictive_search*.
Default identity threshold: 50. (*The user also has the possibility to change the threshold*)  

**Note** - This task can only be performed after you have run OrthoFinder to generate the clusters of orthologs.  



Run ```relaxed_search -h``` to see the usage of this command.

```bash
usage: relaxed_search [-h] -wd WORKINGDIRECTORY -of ORTHOFINDER [-ident IDENT]
                      [-t N_CPU] [-del] [-l] [-v]

Runs the relaxed search step - all the sequences from the database against
representative sequences from the Orthofinder Orthogroups.

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKINGDIRECTORY, --workingDir WORKINGDIRECTORY
                        Working Directory
  -of ORTHOFINDER, --orthofinder ORTHOFINDER
                        Path to OrthoFinder results directory. Inside that
                        directory should be the Orthgroups,
                        Orthogroup_Sequences and WorkingDirectory folders
  -ident IDENT, --identity IDENT
                        Identity threshold to filter the diamond results.
                        DEFAULT: 50
  -t N_CPU              Number of threads to use in the parallel processing.
                        By default it uses all the cpu available on the
                        machine
  -del, --delete        To delete the results stored from diamond (use this
                        option if you don't want to spend memory space,
                        between steps)
  -l, --logfile         To send log messages to a file in the output directory
  -v, --verbose         set loglevel to DEBUG
```

## Restrictive Search  

The second annotation (*restrictive search*) uses as inputs the outputs generated from the *relaxed search*.

Run ```restrictive_search -h``` to see the usage of this command.

```bash
usage: restrictive_search [-h] -wd WORKINGDIRECTORY [-t N_CPU] [-ident IDENT]
                          [-l] [-v]

Runs the restrictive search step - all the sequences from the selected
Orthgroups versus the functions associated with the matched ORAdb hits.

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKINGDIRECTORY, --workingDir WORKINGDIRECTORY
                        Working Directory
  -t N_CPU              Number of threads to use in the parallel processing.
                        By default it uses all the cpu available on the
                        machine.
  -ident IDENT, --identity IDENT
                        Identity threshold to filter the diamond results.
  -l, --logfile         To send log messages to a file in the output directory
  -v, --verbose         set loglevel to DEBUG
```

## Annotation

Annotation of sequences in the clusters of orthologs takes into consideration the following parameters:


**% Identity** - the percentage of identical matches in the range of alignment. Default: 95.

**% Positive Matches** - percent of identical matches + positive matches in the range of alignment.Default: 99.

**% Query Coverage** - percent of the query sequence involved in the range of alignment.Default: 90.

**% Target Coverage** - percent of the target sequence (sequence in the database) involved in the range of alignment.Default: 90.



Run ```annotation -h``` to see the usage of this command.

```bash
usage: annotation [-h] -wd WORKINGDIRECTORY [-s SCORE] [-ident IDENT]
                  [-qc QUERY_COV] [-sc SUBJECT_COV] [-ppos PPOS] [-l] [-v]

Runs the annotation step - filter results from restrictive search and return
the annotated sequences

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKINGDIRECTORY, --workingDir WORKINGDIRECTORY
                        Working Directory
  -s SCORE, --score SCORE
                        score threshold to filter the diamond results.
                        Default: 90
  -ident IDENT, --identity IDENT
                        Identity threshold to filter the diamond results.
  -qc QUERY_COV, --queryCoverage QUERY_COV
                        Query sequence coverage threshold to filter the
                        diamond results.
  -sc SUBJECT_COV, --subjectCoverage SUBJECT_COV
                        Subject sequence coverage threshold to filter the
                        diamond results.
  -ppos PPOS, --percPosMatches PPOS
                        Percentage of positive matches threshold to filter the
                        diamond results (should be higher than identity
                        threshould).
  -l, --logfile         To send log messages to a file int the output
                        directory
  -v, --verbose         set loglevel to DEBUG
```


## Create/Update ORA database (optional)

With this task, you have the option to update your initial database adding the newly annotated sequences from the input genomes.
You have the option to create a new database in a folder of your choice, or simply add the new sequences to the given database.

Run ```create_db -h``` to see the usage of this command.

usage: create_db [-h] -wd WORKINGDIRECTORY (-o OUTPUT | -up) [-l] [-v]

```bash
Runs the create database step - the annotated sequences are added to the
previous database.

optional arguments:
  -h, --help            show this help message and exit
  -wd WORKINGDIRECTORY, --workingDir WORKINGDIRECTORY
                        Working Directory
  -o OUTPUT, --out OUTPUT
                        Output directory to create new database (initial
                        database + new annotated sequences).
  -up, --update         Use this option to update the initial database with
                        the new annotated sequences. ATTENTION: If you use
                        this option the initial database will be changed
                        permanently
  -l, --logfile         To send log messages to a file int the output
                        directory
  -v, --verbose         set loglevel to DEBUG

```


## Running OrtAn

First define the paths of the files and folders you will need:

```bash
work_dir="/path/to/output/folder/" # location where you want to store the results
database="/path/to/database/" # location of the ORAdb (where the FASTA files are located)
orthof="/path/to/orthofinder/results/folder" # location of the folders with the results from OrthoFinder
new_db="/path/to/output/folder/new_db/" # location where the new database should stored (optional)
```

Example:
====

```bash
work_dir="examples/OrtAn_Results/"

database="examples/test_database/"

orthof="examples/" (all folders and files generated from OrthoFinder are located in this folder)
```
Create the necessary directories to store the results:

```bash
mkdir $work_dir
mkdir $new_db (optional)
```

Create a project using as input the database whose path is indicated in the variable database.

All the results and outputs from this project will be stored in the folder indicated in the variable ```$work_dir```.

```bash
create_project -out $work_dir -db $database
```

Run the relaxed_search step, using only 2 cores to run, an identity cutoff of 50% and using as input the OrthoFinder results present in the folder indicated in the variable ```$orthof```.

50% identity cutoff means that all the pair of sequences with an identity percent less than 50 are discarded.

```bash
relaxed_search -wd $work_dir -of $orthof -t 2 -ident 50
```

**Note:** Number of threads should not exceed those available in your machine. If you don't specify the number the tool will use all the available cores in the machine.

Run the restrictive_search step using 2 cores to run.

```bash
restrictive_search -wd $work_dir -t 2
```

Run annotation with thresholds of 95 for identity percent, 99 for positive matches percent and 90 of query and target coverage percent.

```bash
annotation -wd $work_dir -ident 95 -ppos 99 -qc 90 -sc 90 # remove everything after "$work_dir" if default values are to be used
```

Run *create_db* in a new directory adding to the initial database the annotated sequences in the previous step (optional).

```bash
create_db -wd $work_dir -o $new_db
```

## Output

From the *relaxed_search* task, a text file [Associations.txt](examples/OrtAn_Results/Results/Associations.txt) (``` located in examples/OrtAn_Results/Results/Associations.txt```) is generated containing the associations between the clusters of orthologs and the ORAdb functions.

From the *restrictive_search* task we obtain 6 different text files:

[Annotation_Function_Protein.txt](examples/OrtAn_Results/Results/Annotation_Function_Protein.txt) - Shows in the first column the functions and in the second the sequences annotated with those functions (one association per line) (```located in examples/OrtAn_Results/Results/Annotation_Function_Protein.txt```).  

[Annotation_Protein_Function.txt](examples/OrtAn_Results/Results/Annotation_Protein_Function.txt) - Shows in the first column the sequences and in the second the functions assigned (one association per line) (```located in examples/OrtAn_Results/Results/Annotation_Protein_Function.txt```).

[ConOG.txt](examples/OrtAn_Results/Results/ConOG.txt) - Consistent Orthogroups (Clusters of orthologs where all the sequences were annotated to the same function). The function is also indicated (```located in examples/OrtAn_Results/Results/ConOG.txt```).

[DivOG.txt](examples/OrtAn_Results/Results/DivOG.txt) - Divergent Orthogroups (Clusters of orthologs where not all the sequences were annotated to the same function). This means that the ortholog cluster could have sequences that were not annotated to any function or sequences annotated to different functions. These functions are also indicated in the file (```located in examples/OrtAn_Results/Results/DivOG.txt```).

[Orthogroups_Annotation.csv](examples/OrtAn_Results/Results/Orthogroups_Annotation.csv) - This file shows how many sequences in each cluster of orthologs were annotated and to which function (```located in examples/OrtAn_Results/Results/Orthogroups_Annotation.csv```).

[Species_Annotation.csv](examples/OrtAn_Results/Results/Species_Annotation.csv) - This file shows which functions are present in which species (1 - at least one sequence of a species annotated to the function, 0 - no sequences annotated to the function) (```located in examples/OrtAn_Results/Results/Species_Annotation.csv```).


Step 5) Interspecies interactions
==== 

Identification of interspecies interactions can be performed only taking into account the individual species' functional potential or by adding constraints related to reaction subsets that need to be performed by single species or even the presence of transport reactions.

Requirements include:

-> [final_gpr.xlsx](examples/test_database/final_gpr.xlsx) - generated from running the *get_gpr.sh* script

-> [Species_Annotation.csv](microbial_interactions/data/Species_Annotation.csv) - generated from OrtAn

-> [user_input.csv](examples/OrtAn_Results/Results/test_user_input.csv) - a file provided by the user defining the constraints of each pathway or sets of reactions of interest contained in ORAdb 

-> [GP_rules.json](microbial_interactions/data/GP_rules.json) - generated from the tool

-> [paths.json](microbial_interactions/data/paths.json) - generated from the tool 

-> [species_to_exclude](microbial_interactions/data/species_to_exclude.json) - generated from the tool

-> module_list.txt / pathway_list.txt - Additional subsetting of reactions from the ORAdb (optional)


Identification of species with the genetic potential to perform complete pathways individually or by interactions with other requires the execution of the commands shown below.


**Command 1: Depending on the subsetting of reactions from ORAdb you can use the following (choose only one option):**


*Subsetting to a pathway* - takes as inputs the pathway_list, final_grp.xlsx, Species_Annotation.csv and user_input.csv files. Output_folder is defined by the user.

> Rscript /path/to/folder/gpr_manipulation.R  -p /path/to/folder/pathway_list.txt -n /path/to/folder/final_gpr.xlsx -s /path/to/folder/Species_Annotation.csv -u /path/to/folder/user_input.csv -o /path/to/output_folder

*Subsetting to a module list*

> Rscript /path/to/folder/gpr_manipulation.R -m /path/to/folder/module_list.txt -n /path/to/folder/final_gpr.xlsx -s /path/to/folder/Species_Annotation.csv -u /path/to/folder/user_input.csv -o /path/to/output_folder

*Using the complete ORAdb (DEFAULT)*

>Rscript /path/to/folder/gpr_manipulation.R -n /path/to/folder/final_gpr.xlsx -s /path/to/folder/Species_Annotation.csv -u /path/to/folder/user_input.csv -o /path/to/output_folder

During this task the *GP_rules.json*, *paths.json* and *species_to_exclude.json* are also generated by the script and stored in the output_folder defined by the *-o* flag.

**Command 2: Adding constraints to reduce search space:**

```
sh combinations.sh /path/to/folder/Species_Annotation.csv /path/to/folder/GP_rules.json /path/to/folder/paths.json /path/to/folder/species_to_exclude.json > /path/to/folder/output_file_combinatinations.txt
```

where *output_file_combinations.txt* is a name defined by the user.

Citing OrtSuite
====

OrtSuite will be submitted to BioxRiv by the end of July and a link will be added at that time. If other software contained and used by OrtSuite was also useful in your research (e.g. DIAMOND, BLAST and OrthoFinder) please give them credit as well.

Contributions
====

Authors of pipeline: Joao Saraiva and Ulisses Nunes da Rocha.

Principal Investigator: Ulisses Nunes da Rocha

Institution: Microbial Data Sciences group, Helmholtz Center for Environmental Research, Department of Environmental Microbiology, Leipzig, Germany

All feedback is welcome. For errors and bugs, please open a new Issue thread on this github page, and we will try to address them as soon as possible. For general feedback you can contact us at mds@ufz.de.

