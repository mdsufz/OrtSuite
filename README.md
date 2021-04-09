# OrtSuite - a flexible pipeline for annotation of ecosystem processes and prediction of putative microbial interactions

OrtSuite was developed with the goal to facilitate annotation of ecosystem processes and identify putative microbial interactions by automating the complete process from sequence retrieval, clustering of ortholog sequences, functional annotation, to putative microbial interactions prediction. 



![workflow](https://github.com/mdsufz/OrtSuite/blob/master/workflow_ortSuite.png)

**OrtSuite workflow** 

**(a)** Protein sequence from samples supplied by the user are clustered using OrthoFinder. 

**(b)** A text file is provided containing a list of identifiers for each reaction in the pathway of interest supplied by the user to retrieve all protein sequences from KEGG.

**(c)** Sequences mapped to reactions are stored in ORAdb.

**(d)** Functional annotation consists of a two-stage process (relaxed and restrictive search). Relaxed search **(e)** performs sequence alignments between 50% of randomly selected sequences from each generated cluster. Clusters whose representative sequences share a minimum 0.001 E-value to sequences in reaction set(s) in ORAdb transition to the restrictive search **(f)**. Here, all sequences from the cluster is aligned to all sequences in the corresponding reaction set(s) to which they had a hit (Default E-value threshold of 1E-9). 

**(g and h)** The annotated sequences are used to identify putative microbial interactions based on their functional potential.

**(rounded blue rectangles)** Additional constraints can be added to reduce the search space of microbial interactions.


# Overview of OrtSuite

**Installation**

**ORAdb construction and GPR definition:** Generation of the user-defined Ortholog-Reaction Association (ORAdb) database and download of Gene-Protein-Reaction rules from KEGG.

**Clustering of orthologs**   

**Functional annotation of clusters of orthologs**   

**Prediction of interspecies interactions based on the functional potential of individual species**   


# System Requirements

Resources for Ortsuite will vary depending on the amount of data being processed. In the example provided (consisting of 7 reactions with 16 associated KEGG Ortholog identifiers), we used an Intel Core i5-6200U 2.3GHz with 4 cores and 16 Gb of RAM. OrtSuite officially supports only Linux OS. 

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Installation
============


# OrtSuite 

OrtSuite is a python tool that performs functional annotation of clusters of orthologs and identifies putative microbial interactions. This tool automatically retrieves sequence data in bulk from KEGG (Kyoto Encyclopedia of Genes and Genomes) database to generate the Ortholog Reaction-Associated user-defined database (*ORAdb*).
Generation of clusters of orthologs is performed by OrthoFinder.

**Requirements:**  Python 3.6

**Dependencies:**  Setuptools, bs4, grequests, [OrthoFinder](https://github.com/davidemms/OrthoFinder), [DIAMOND](https://github.com/bbuchfink/diamond)


We suggest the use of conda to install the virtual environment.

Create a virtual environment.

```bash
conda create -n OrtSuite python=3.6

```
Activate the virtual environment.

```bash
conda activate OrtSuite
```

Install dependencies

```bash
pip install -U Setuptools
pip install grequests
pip install BeautifulSoup4=4.8.2
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



## Steps 1 and 2) Generating the Ortholog Reaction-Assosiation database (ORAdb) and download of Gene-Protein-Reaction (GPR) rules from KEGG.


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

ORAdb Output
======
A FASTA format file for each selected KO will be placed in the output folder defined by the user.
If you use -e or -r option, an additional file will be generated *associations.txt*, which indicates which kos where selected for download for each reaction/EC number.
In the same folder a file *info_db.csv* contains a table with information regarding the KO's that were selected for download along with their name and associated EC numbers and Reactions IDs. 
In the case of using KO identifiers the file associations.txt must be manually added (for an example please see [associations.txt](examples/associations.txt)).

GPR rules
=====
Once the FASTA files containing the sequences for the list of KOs is completed the GPR rules are obtained.

**Note:** The gpr rules associated with each identifier is stored in a file called *final_gpr.xlsx*. 


**Both the ORAdb construction and download of the GPR rules is performed by running the script *DB_construction.sh*.**

This script takes as inputs:

1 - project full path (It needs to be inside the OrtSuite main folder. Ex: "~/OrtSuite/PROJECT_NAME/")
2 - the flag of the list (-m, -r, -e, or -k) 
3 - full path of the list of identifiers

Examples of commands: 
====
When using a single pathway map from KEGG
```bash
sh DB_construction.sh ~/OrtSuite/PROJECT_NAME/ -m /path/to/folder/map00362
```
When using reaction identifiers
```bash
sh DB_construction.sh ~/OrtSuite/PROJECT_NAME/ -r /path/to/folder/reaction_list.txt
```
When using EC numbers
```bash
sh DB_construction.sh ~/OrtSuite/PROJECT_NAME/ -e /path/to/folder/ec_list.txt
```
When using KO identifiers
```bash
sh DB_construction.sh ~/OrtSuite/PROJECT_NAME/ -k /path/to/folder/ko_list.txt

```
**An output_folder called *database* will be created in the project directory that the user provided.**

To test if the tool is working you can use the files contained in the [examples](https://github.com/mdsufz/OrtSuite/blob/master/examples) folder.

Note: Running this command may take some time and memory space.

**Due to the limited information in KEGG database concerning GPR rules, manual inspection of the *final_gpr.xlsx* is strongly recommended.
The script generates a table with the GPR rules for all reaction-enzyme pairs. Since the same reaction can occur in different modules with different gene rules, the user should edit this file so that it only includes ONE unique rule per REACTION according to the target pathway.
For a more comprehensive explanation please see the [tutorial](OrtSuite_tutorial.md) with an example.**  

Step 3) Clustering of orthologs (OrthoFinder)
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

Steps 4 and 5) Functional annotation of clusters of orthologs and prediction of interspecies interactions
====

Functional annotation of clusters of orthologs generated from **Step 3** is based on the user-defined ORAdb generated in **Step 1**.

Overview of functional annotation:

>Create Project: receive the information of the ORAdb database
>               and create the working directory structure necessary for the subsequent tasks.

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



**Functional annotation and prediction of interspecies interactions is performed by running the script *annotate_and_predict.sh*.**

This script takes as inputs:

1 - project full path (It needs to be inside the OrtSuite main folder. Ex: "~/OrtSuite/PROJECT_NAME/")
2 - full path to the results from Orthofinder (Step 3) (Ex: "/path/to/output/folder/Results")
3 - full path of the text file with the user-defined constraints (Ex: -> [user_input.csv](examples/OrtAn_Results/Results/test_user_input.csv))



**Before running OrtAn you need to:**

- Run OrthoFinder with the input genomes  
- Prepare the database with the necessary format.  


Example:
====

```bash
sh annotate_and_predict.sh ~/OrtSuite/PROJECT_NAME/ ~/path/to/folder/orthofinder/Results_Jan31/ /path/to/folder/user_input.csv
```

## Functional annotation output files

From the *relaxed_search* task, a text file [Associations.txt](examples/OrtAn_Results/Results/Associations.txt) (``` located in examples/OrtAn_Results/Results/Associations.txt```) is generated containing the associations between the clusters of orthologs and the ORAdb functions.

From the *restrictive_search* task we obtain 6 different text files:

[Annotation_Function_Protein.txt](examples/OrtAn_Results/Results/Annotation_Function_Protein.txt) - Shows in the first column the functions and in the second the sequences annotated with those functions (one association per line) (```located in examples/OrtAn_Results/Results/Annotation_Function_Protein.txt```).  

[Annotation_Protein_Function.txt](examples/OrtAn_Results/Results/Annotation_Protein_Function.txt) - Shows in the first column the sequences and in the second the functions assigned (one association per line) (```located in examples/OrtAn_Results/Results/Annotation_Protein_Function.txt```).

[ConOG.txt](examples/OrtAn_Results/Results/ConOG.txt) - Consistent Orthogroups (Clusters of orthologs where all the sequences were annotated to the same function). The function is also indicated (```located in examples/OrtAn_Results/Results/ConOG.txt```).

[DivOG.txt](examples/OrtAn_Results/Results/DivOG.txt) - Divergent Orthogroups (Clusters of orthologs where not all the sequences were annotated to the same function). This means that the ortholog cluster could have sequences that were not annotated to any function or sequences annotated to different functions. These functions are also indicated in the file (```located in examples/OrtAn_Results/Results/DivOG.txt```).

[Orthogroups_Annotation.csv](examples/OrtAn_Results/Results/Orthogroups_Annotation.csv) - This file shows how many sequences in each cluster of orthologs were annotated and to which function (```located in examples/OrtAn_Results/Results/Orthogroups_Annotation.csv```).

[Species_Annotation.csv](examples/OrtAn_Results/Results/Species_Annotation.csv) - This file shows which functions are present in which species (1 - at least one sequence of a species annotated to the function, 0 - no sequences annotated to the function) (```located in examples/OrtAn_Results/Results/Species_Annotation.csv```).

## Interspecies interactions output files

Interspecies interactions will be stored in a folder called *interactions* located inside the project folder (Ex: /project_folder/work_dir/interactions).

The complete list of files generated during prediction of interspecies interactions are:


- [complete_pathway_species.txt](complete_pathway_species.txt) : where all species with the complete functional potential for each pathway of interest are listed.
- [Reactions_mapped_to_species.csv](Reactions_mapped_to_species.csv) : A binary table showing species that possess the genomic content to encode proteins involved in each reaction in ORAdb (1- present, 0 - absent).
- [single_org_subset_interactions.txt](single_org_subset_interactions.txt) :  where the interspecies interactions that fulfill the constraint of reaction subsets required to be performed by individual species (e.g. reactions X and Y have to be present in a single organism) are shown.
- Files containing the number of species with the functional potential to each reaction. For example, [Aerobic conversion of benzoate to acetyl-CoA_species_per_reactions.txt](Aerobic conversion of benzoate to acetyl-CoA_species_per_reactions.txt)
- A file containing all interspecies interactions whose combined functional potential allow a complete pathway of interest (defined in *user_input.csv*). As an example: [Aerobic conversion of benzoate to acetyl-CoA](Aerobic_benzoate-acetylCoA.csv)


Citing OrtSuite
====

A preprint version of OrtSuite is available [here](https://www.researchsquare.com/article/rs-52281/v1). A link to the published manuscript will be provided as soon as possible. If other software contained and used by OrtSuite was also useful in your research (e.g. DIAMOND, BLAST and OrthoFinder) please give them credit as well.

Contributions
====

Authors of pipeline: Joao Saraiva and Ulisses Nunes da Rocha.

Principal Investigator: Ulisses Nunes da Rocha

Institution: Microbial Data Sciences group, Helmholtz Center for Environmental Research, Department of Environmental Microbiology, Leipzig, Germany

All feedback is welcome. For errors and bugs, please open a new Issue thread on this github page, and we will try to address them as soon as possible. For general feedback you can contact us at mds@ufz.de.

