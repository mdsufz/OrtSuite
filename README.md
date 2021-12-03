# OrtSuite - a flexible pipeline for annotation of ecosystem processes and prediction of putative microbial interactions

OrtSuite was developed with the goal to facilitate annotation of ecosystem processes and identify putative microbial interactions by automating the complete process from sequence retrieval, clustering of ortholog sequences, functional annotation, to putative microbial interactions prediction. OrtSuite only requires three commands to run (from database construction to clustering of orthologs, functional annotation and putative microbial synergistic interactions). 
We provide video tutorials briefling explaining how to use OrtSuite [here](https://www.youtube.com/playlist?list=PLH4_1OSQOdnqTQY-olYX3vZss4DWtPkZj).



![workflow](https://github.com/mdsufz/OrtSuite/blob/master/Figure_1.png)

**OrtSuite workflow** 

**(Task 1)** OrtSuite takes a text file containing a list of identifiers for each reaction in the pathway of interest supplied by the user to retrieve all protein sequences from KEGG Orthology and are stored in ORAdb. Subsequently the same list of identifiers is used to obtain the Gene-Protein-Reaction (GPR) rules from KEGG Modules.

**(Task 2)** Protein sequences, supplied by the user, are clustered using OrthoFinder.

**(Task 3)** Functional annotation, identification of synergistic species interactions and generation of a graphical visualization of the network. Functional annotation consists of a two-stage process (relaxed and restrictive search). Relaxed search performs sequence alignments between 50% of randomly selected sequences from each generated cluster. Clusters whose representative sequences share a minimum E-value of 0.001 to sequences in the reaction set(s) in ORAdb transition to the restrictive search . Here, all sequences from the cluster are aligned to all sequences in the corresponding reaction set(s) to which they had a hit (default E-value = 1e-9). Annotated sequences are further filtered to those with a bit score greater than 50. The annotated sequences are used to identify putative microbial interactions based on their functional potential. Additional constraints can be added to reduce the search space of microbial interactions (e.g. subsets of reactions required to be performed by single species, transport-related reactions). A graphical network based on the reaction set defined by user allows to visualize and filter, interactively, the results obtained.


# Overview of OrtSuite

**Installation**

**Task 1) ORAdb construction and GPR definition:** Generation of the user-defined Ortholog-Reaction Association (ORAdb) database and download of Gene-Protein-Reaction rules from KEGG.

**Task 2) Clustering of orthologs**   

**Task 3) Functional annotation of clusters of orthologs and Prediction of interspecies interactions based on the functional potential of individual species**   


# System Requirements

Resources for Ortsuite will vary depending on the amount of data being processed. In the example provided (consisting of 7 reactions with 16 associated KEGG Ortholog identifiers), we used an Intel Core i5-6200U 2.3GHz with 4 cores and 16 Gb of RAM. OrtSuite officially supports only Linux OS. 

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Installation
============


# OrtSuite 

OrtSuite is a python tool that performs functional annotation of clusters of orthologs and identifies putative microbial interactions. This tool automatically retrieves sequence data in bulk from KEGG (Kyoto Encyclopedia of Genes and Genomes) database to generate the Ortholog Reaction-Associated user-defined database (*ORAdb*).
Generation of clusters of orthologs is performed by OrthoFinder.

With Docker image
=====

For personal computers or High Perfomance Computers (HPCs) that accept docker images we recommend the following intallation guide:


**Requirements:** [Docker](https://docs.docker.com/engine/install/ubuntu/)  
**To note: Installing via docker does not require you to install the remaining dependencies (e.g. OrthoFinder, pandoc)!!**

Run the following command to pull the docker image

```bash
sudo docker pull mdsufz/ortsuite:latest
```
To execute the docker image run the following command in your terminal

```bash
sudo docker run -it --name ortsuite_docker mdsufz/ortsuite bash
cd OrtSuite
```
***Note***: The location of OrthoFinder using the docker image is *app/OrthoFinder*.

With Conda install
=====
**If you do not wish to use a docker image the following procedure is required.**

**Requirements:**  Python 3.6

**Dependencies:**  [OrthoFinder](https://github.com/davidemms/OrthoFinder), [DIAMOND](https://github.com/bbuchfink/diamond), R, [pandoc](https://pandoc.org/installing.html)  

**First you need to clone the OrtSuite [repository](https://github.com/mdsufz/OrtSuite)**

We suggest the use of miniconda to install the virtual environment.

Get Miniconda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Run the following to allow for execution of Miniconda

```bash
chmod +x Miniconda3-latest-Linux-x86_64.sh
sh ./Miniconda3-latest-Linux-x86_64.sh
export PATH=~/miniconda/bin:$PATH
```


Create a virtual environment.

```bash
conda env create -f ortsuite_env.yml -n ortsuite_env

```
Activate the virtual environment.

```bash
conda activate ortsuite_env
```

Next move to the folder where the file setup.py from OrtSuite is located.

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


Pandoc
====

This library is needed to generate the interactive network.
```bash
sudo apt-get install pandoc
```

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Using OrtSuite
=====


Once installation of OrtSuite and all dependencies are completed the different commands can be called independently.



## Task 1) Generating the Ortholog Reaction-Assosiation database (ORAdb) and download of Gene-Protein-Reaction (GPR) rules from KEGG.


The generation of a user-defined ortholog reaction-association database starts with the download of all the sequences associated with the KO (KEGG Orthology) group to a FASTA file.


The input can be:

- a KEGG pathway map ID (e.g. map00362)

- [List of KO IDs](https://github.com/mdsufz/OrtSuite/blob/master/examples/kos.txt)

- [List of KEGG Reaction IDs](https://github.com/mdsufz/OrtSuite/blob/master/examples/rx.txt)

- [List of EC (Enzyme commission) numbers](https://github.com/mdsufz/OrtSuite/blob/master/examples/ecs.txt)

Note: The format of the input lists must be a txt file with only one ID per line. 


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

1 - project full path ("path/to/folder/PROJECT_NAME/")  
2 - the flag of the list (-m, -r, -e, or -k)  
3 - full path of the list of identifiers  
4 - full path of the OrtSuite installation firectory (e.g. "~/OrtSuite/")  

Examples of commands: 
====
When using a single pathway map from KEGG
```bash
sh DB_construction.sh /PROJECT_NAME/ -m /path/to/folder/map00362 ~/OrtSuite
```
When using reaction identifiers
```bash
sh DB_construction.sh /PROJECT_NAME/ -r /path/to/folder/reaction_list.txt ~/OrtSuite
```
When using EC numbers
```bash
sh DB_construction.sh /PROJECT_NAME/ -e /path/to/folder/ec_list.txt ~/OrtSuite
```
When using KO identifiers
```bash
sh DB_construction.sh /PROJECT_NAME/ -k /path/to/folder/ko_list.txt ~/OrtSuite

```
**An output_folder called *database* will be created in the project directory that the user provided.**

To test if the tool is working you can use the files contained in the [examples](https://github.com/mdsufz/OrtSuite/blob/master/examples) folder.

Note: Running this command may take some time and memory space.

**Due to the limited information in KEGG database concerning GPR rules, manual inspection of the *final_gpr.xlsx* is strongly recommended.
The script generates a table with the GPR rules for all reaction-enzyme pairs. Since the same reaction can occur in different modules with different gene rules, the user should edit this file so that it only includes ONE unique rule per REACTION according to the target pathway.


Task 2) Clustering of orthologs (OrthoFinder)
====

OrthoFinder takes as input a folder containing the FASTA sequences the user wants to cluster.

```bash
~/OrthoFinder/orthofinder -f ~/path/to/sequence/folder -o /path/to/output/folder -og
```

Note: If you wish to use BLAST+ instead of DIAMOND please use the following:

```bash
~/OrthoFinder/orthofinder -f ~/path/to/sequence/folder -o /path/to/output/folder -S blast -og
```

Note: If you wish to use MMSeqs2 instead of DIAMOND please use the following:

```bash
~/OrthoFinder/orthofinder -f ~/path/to/sequence/folder -o /path/to/output/folder -S mmseqs -og
```

**Note:** OrthoFinder's output folder is generated automatically. The user can, however, define the parent directory where to store the output folder (e.g. /Documents/).

Task 3) Functional annotation of clusters of orthologs, prediction of interspecies interactions and graphical visualization of the network
====

Functional annotation of clusters of orthologs generated from **Task 2** is based on the user-defined ORAdb generated in **Task 1**.

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



**All steps included in Task 3 are performed by running the script *annotate_and_predict.sh*.**

This script takes as inputs:

1 - project full path ("path/to/folder/PROJECT_NAME/")  
2 - full path to the results from Orthofinder (Step 3) ("/path/to/output/folder/Results")  
3 - full path of the text file with the user-defined constraints (Ex: [user_input.csv](examples/OrtAn_Results/Results/test_user_input.csv))  
4 - full path of the OrtSuite installation firectory (e.g. "/OrtSuite/")  
5 - full path of the list of reaction pairs (separated by tab) for a directed network (optional) (e.g. "/OrtSuite/reaction_pairs.tsv")  
  

**Before running this script you need to:**

- Run OrthoFinder with the input genomes  
- Prepare the database with the necessary format.  


Example:
====

```bash
sh annotate_and_predict.sh /PROJECT_NAME/ ~/path/to/folder/orthofinder/Results_Jan31/ /path/to/folder/user_input.csv ~/OrtSuite/
```

## Functional annotation output files

From the *relaxed_search* task, a text file [Associations.txt](examples/OrtAn_Results/Results/Associations.txt) (``` located in examples/OrtAn_Results/Results/Associations.txt```) is generated containing the associations between the clusters of orthologs and the ORAdb functions.

From the *restrictive_search* task we obtain 6 different text files:

- [Annotation_Function_Protein.txt](examples/OrtAn_Results/Results/Annotation_Function_Protein.txt) - Shows in the first column the functions and in the second the sequences annotated with those functions (one association per line) (```located in examples/OrtAn_Results/Results/Annotation_Function_Protein.txt```).  

- [Annotation_Protein_Function.txt](examples/OrtAn_Results/Results/Annotation_Protein_Function.txt) - Shows in the first column the sequences and in the second the functions assigned (one association per line) (```located in examples/OrtAn_Results/Results/Annotation_Protein_Function.txt```).

- [ConOG.txt](examples/OrtAn_Results/Results/ConOG.txt) - Consistent Orthogroups (Clusters of orthologs where all the sequences were annotated to the same function). The function is also indicated (```located in examples/OrtAn_Results/Results/ConOG.txt```).

- [DivOG.txt](examples/OrtAn_Results/Results/DivOG.txt) - Divergent Orthogroups (Clusters of orthologs where not all the sequences were annotated to the same function). This means that the ortholog cluster could have sequences that were not annotated to any function or sequences annotated to different functions. These functions are also indicated in the file (```located in examples/OrtAn_Results/Results/DivOG.txt```).

- [Orthogroups_Annotation.csv](examples/OrtAn_Results/Results/Orthogroups_Annotation.csv) - This file shows how many sequences in each cluster of orthologs were annotated and to which function (```located in examples/OrtAn_Results/Results/Orthogroups_Annotation.csv```).

- [Species_Annotation.csv](examples/OrtAn_Results/Results/Species_Annotation.csv) - This file shows which functions are present in which species (1 - at least one sequence of a species annotated to the function, 0 - no sequences annotated to the function) (```located in examples/OrtAn_Results/Results/Species_Annotation.csv```).

## Interspecies interactions output files and network visualization

Interspecies interactions will be stored in a folder called *interactions* located inside the project folder (Ex: /project_folder/work_dir/interactions).

The complete list of files generated during prediction of interspecies interactions are:


- [complete_pathway_species.txt](complete_pathway_species.txt) : where all species with the complete functional potential for each pathway of interest are listed.
- [Reactions_mapped_to_species.csv](Reactions_mapped_to_species.csv) : A binary table showing species that possess the genomic content to encode proteins involved in each reaction in ORAdb (1- present, 0 - absent).
- [single_org_subset_interactions.txt](single_org_subset_interactions.txt) :  where the interspecies interactions that fulfill the constraint of reaction subsets required to be performed by individual species (e.g. reactions X and Y have to be present in a single organism) are shown.
- Files containing the number of species with the functional potential to each reaction. For example, [Aerobic conversion of benzoate to acetyl-CoA_species_per_reactions.txt](Aerobic conversion of benzoate to acetyl-CoA_species_per_reactions.txt)
- A file containing all interspecies interactions whose combined functional potential allow a complete pathway of interest (defined in *user_input.csv*). As an example: [Aerobic conversion of benzoate to acetyl-CoA](Aerobic_benzoate-acetylCoA.csv)
- A HTML file containing the interactive network visulization for the pathway of interest. For example: [network_example](network_example.png). 

Citing OrtSuite
====

Publication of OrtSuite is available [here](https://www.life-science-alliance.org/content/4/12/e202101167). If other software contained and used by OrtSuite was also useful in your research (e.g. DIAMOND, BLAST and OrthoFinder) please give them credit as well.

Contributions
====

Authors of pipeline: Joao Saraiva and Ulisses Nunes da Rocha.

Principal Investigator: Ulisses Nunes da Rocha

Institution: Microbial Data Sciences group, Helmholtz Center for Environmental Research, Department of Environmental Microbiology, Leipzig, Germany

All feedback is welcome. For errors and bugs, please open a new Issue thread on this github page, and we will try to address them as soon as possible. For general feedback you can contact us at mds@ufz.de.

