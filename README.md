# OrtSuite - a flexible pipeline for annotation of ecosystem processes and prediction of putative microbial interactions

OrtSuite was developed with the goal to facilitate annotation of ecosystem processes and identify putative microbial interactions by automating the complete process from sequence retrieval, clustering of ortholog sequences, functional annotation, to putative microbial interactions prediction. 



![workflow](https://github.com/msdsufz/OrtSuite/blob/master/workflow_ortan_no_legend.png)
**OrtAn workflow** - Protein sequence from samples supplied by the user are clustered using OrthoFinder **(a)**. OrtScraper takes a text file containing a list of identifiers for each reaction in the pathway of interest supplied by the user to retrieve all protein sequences from KEGG **(b)**. Sequences mapped to reactions are stored in ORAdb **(c)**. Functional annotation **(d)** consists of a two-stage process (relaxed and restrictive search). Relaxed search **(e)** performs sequence alignments between 10% of randomly selected sequences from each generated cluster. Clusters whose representative sequences share a minimum 50% identity to sequences in reaction set(s) in ORAdb transition to the restrictive search **(f)**. Here, all sequences from the cluster is aligned to all sequences in the corresponding reaction set(s) to which they had a hit. Finally, the annotated sequences are used to identify putative microbial interactions based on their functional potential **(g and h)**. Additional constraints can be added to reduce the search space of microbial interactions **(rounded blue rectangles)**.


# Overview of OrtSuite


**OrthoFinder:** Clustering of orthologs

**OrtScraper:** Bulk download of protein sequences for populating a user-defined database

**OrtAn:** Functional annotation of clusters of orthologs


# System Requirements

Resources for Ortsuite will vary depending on the amount of data being processed. In the example provided (consisting of 7 reactions with 16 associated KEGG Ortholog identifiers), we used an Intel Core i5-6200U 2.3GHz with 4 cores and 16 Gb of RAM. OrtSuite officially supports only Linux OS. 



Installation
============


# OrtSuite

OrtSuite is a python tool that performs functional annotation of clusters of orthologs and identifies putative microbial interactions. This tool automatically retrieves sequence data in bulk from KEGG (Kyoto Encyclopedia of Genes and Genomes) database to generate the Ortholog Reaction-Associated user-defined database (*ORAdb*).
Generation of clusters of orthologs is performed by OrthoFinder.

**Requirements:**  Python 3

**Dependencies:**  bs4, grequests, [OrthoFinder](https://github.com/davidemms/OrthoFinder), [DIAMOND](https://github.com/bbuchfink/diamond)


We suggest the use of the package manager [pip](https://pip.pypa.io/en/stable/) to install the virtual environment (*virtualenv*).


```bash
pip install virtualenv
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
python3 setup.py install
```

Once the installation is finished the tool should be ready to use.


OrthoFinder
=====

Install OrthoFinder (Available here: https://github.com/davidemms/OrthoFinder)
```bash
Download latest version: wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz 
Decompress file: tar xzf OrthoFinder.tar.gz (to your desired directory)
Testing: ./OrthoFinder/orthofinder -h (help menu)
```
Dependencies of OrthoFinder

*MCL*

The mcl clustering algorithm is available in the repositories and can be installed as any other package.
```bash
sudo apt-get install mcl
```

*Note* - No other dependency of OrthoFinder is required to run OrtSuite since we only use it to perform the clustering of orthologs. 


DIAMOND
====

Install DIAMOND (Available here: https://github.com/bbuchfink/diamond/releases)
```bash
Download latest version: wget https://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz
Decompress file: tar xzf diamond-linux64.tar.gz
Copy DIAMOND to your local directory: sudo cp diamond /usr/local/bin
```



Usage
=====

Please view the OrtSuite [tutorial](https://github.com/msdsufz/OrtSuite/OrtSuite_tutorial.md) for detailed instructions and examples.

Once installation of OrtSuite and all dependencies are completed the different commands can be called independently.



## download_kos

Download all the sequences associated with the KO (KEGG Orthology) group to a FASTA file.

The input can be:

- a KEGG pathway map ID

- List of KO IDs

- List of KEGG Reaction IDs

- List of EC (Enzyme commission) numbers

Note: The format of the input lists must be a txt file with only one ID per line. In the case of using KO identifiers the file associations.txt must be manually added (for an example please see [associations.txt](https://github.com/msdsufz/OrtSuite/examples/associations.txt))


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

To test if the tool is working you can use the files contained in the examples folder.

Note: Running this command may take some time and memory space.


download_kos Output
======
In the output folder you will find one FASTA file for each selected KO.
If you use -e or -r option, an additional file will be generated *associations.txt*, which indicates which kos where selected for download for each reaction/EC number.
In the same folder a file *info_db.csv* contains a table with information regarding the KO's that were selected for download along with their name and associated EC numbers and Reactions IDs.



OrthoFinder
====

OrthoFinder takes as input a folder containing the FASTA sequences the user wants to cluster.

```bash
~/OrthoFinder/orthofinder -f ~/path/to/sequence/folder -og
```


OrtAn
====

This tool uses the output from [OrthoFinder](https://github.com/davidemms/OrthoFinder) to perform the annotation of the generated clusters of orthologs based on the user-defined ORAdb.

Overview of OrtAn
```bash
Create Project: OrtAn receives the information of the ORAdb database and creates the working directory structure necessary for the subsequent tasks.
Relaxed Search: This task takes the OrthoFinder information and runs a first relaxed search with DIAMOND to identify the associations between the returned clusters of orthologs and the different functions of the ORA database.
Restrictive Search: This task performs a restrictive search only between the clusters of orthologs and groups of functions from the ORA database that were related during the relaxed_search.
Annotation: This task consists in the annotation of the sequences present in the clusters of orthologs after performing the restrictive_search.
Identification of putative microbial interactions: This task allows the user to extract the putative microbial interactions based on different sets of constraints (e.g. number of interacting species, ability to perform complete pathways, groups of interactions where a single species is responsible for a subset of reactions in the pathway).
```

### Inputs:

- OrthoFinder output directories

- ORAdb


# Pipeline

Before running OrtAn you need to:
- Run OrthoFinder with the input genomes;
- Prepare the database with the necessary format.

### OrthoFinder required files and directories:

```/Orthogroups/Orthogroups.txt```

```/WorkingDirectory/SpeciesIDs.txt```

```/WorkingDirectory/SequenceIDs.txt```

```/Orthogroup_Sequences/```


## Create Project


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

**Note** - This task can only be performed after you have run OrthoFinder to generate the clusters of orthologs.

Ten (10) percent of random sequences from each cluster is used to retrieve the associations which will reduce the search space on the second annotation *restrictive_search*.
Default identity threshold: 50. (*The user also has the possibility to change the threshold*) 

Run ```relaxed_search -h``` to see the usage of this command.

```bash
usage: relaxed_search [-h] -wd WORKINGDIRECTORY -of ORTHOFINDER [-ident IDENT]
                      [-t N_CPU] [-del] [-l] [-v]

Runs the relaxed search step - all the sequences from the database against
representative sequences from the Orthofinder Orthgroups.

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


Run ```restrictive_search -h``` to see the usage of this command.

```bash
usage: restrictive_search [-h] -wd WORKINGDIRECTORY [-t N_CPU] [-ident IDENT]
                          [-l] [-v]

Runs the restrictive search step - all the sequences from the selected
Orthgroups versus the functions associated with.

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
```bash

**% Identity** - the percentage of identical matches in the range of alignment. Default: 95.

**% Positive Matches** - percent of identical matches + positive matches in the range of alignment.Default: 99.

**% Query Coverage** - percent of the query sequence involved in the range of alignment.Default: 90.

**% Target Coverage** - percent of the target sequence (sequence in the database) involved in the range of alignment.Default: 90.
```

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

## Create Database

With this command, you have the option to update your initial database adding the newly annotated sequences from the input genomes.
You can have the option to create a new database in a folder of your choice, or simply add the new sequences to the given database.

Run ```create_db -h``` to see the usage of this command.




## Run the pipeline

**Note:** You can try this with the example data.

Define the important paths you need:

```bash
work_dir="/path/to/output/folder/"
database="/path/to/database/"
orthof="/path/to/orthofinder/results/folder"
new_db="/path/to/output/folder/new_db/"
```

Create the necessary directories to store the results:

```bash
mkdir $work_dir
mkdir $new_db
```

Run all the commands in order.

Create a project using as input the database which the path is indicated in the variable database.

All the results and outputs from this project will be stored in the folder indicated in the variable ```$work_dir```.

```bash
create_project -out $work_dir -db $database
```

Run the relaxed_search step, using only 2 cores to run, an identity cutoff of 80% and using as input the OrthoFinder results present in the folder indicated in the variable ```$ortho```.

80% identity cutoff means that all the pair of sequences with an identity percent less than 80 are discarded.

```bash
relaxed_search -wd $work_dir -of $orthof -t 2 -ident 80
```

**Note:** Number of threads should not exceed the one available in your machine. If you don't specify the number the tool will use all the available cores in the machine.

Run the restrictive_search step using 2 cores to run.

```bash
restrictive_search -wd $work_dir -t 2
```

Run annotation with thresholds of 95 for identity percent, 99 for positive matches percent and 90 of query and target coverage percent.

```bash
annotation -wd $work_dir -ident 95 -ppos 99 -qc 90 -tc 90
```

Run create_db in a new directory (the one indicated in the variable new_db) adding to the initial database the annotated sequences in the previous step.

```bash
create_db -wd $work_dir -o $new_db
```

## Output

From relaxed_search step, a text file (```/Results/Associations.txt```) is generated containing the associations between the clusters of orthologs and the ORAdb functions.

From the annotation we obtain 6 different text files:

```/Results/Annotation_Function_Protein.txt``` - Shows in the first column the functions and in the second the sequences annotated with those functions (one association per line).

```/Results/Annotation_Protein_Function.txt``` - Shows in the first column the sequences and in the second the functions assigned (one association per line).

```/Results/ConOG.txt``` - Consistent Orthogroups (Clusters of orthologs where all the sequences were annotated to the same function). The function is also indicated.

```/Results/DivOG.txt``` - Divergent Orthogroups (Clusters of orthologs where not all the sequences were annotated to the same function). This means that the ortholog cluster could have sequences that were not annotated to any function or sequences annotated to different functions. These functions are also indicated in the file.

```/Results/Orthogroups_Annotation.csv``` - This file shows how many sequences in each cluster of orthologs were annotated and to which function.

```/Results/Species_Annotation.csv``` - This file shows which functions are present in which species (1 - at least one sequence of a species annotated to the function, 0 - no sequences annotated to the function).




