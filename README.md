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

Install DIAMOND ((Available here: https://github.com/bbuchfink/diamond/releases)
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

Note: The format of the input lists must be a txt file with only one ID per line. In the case of using KO identifiers the file associations.txt must be manually added (for an example please see [associations](https://github.com/msdsufz/OrtSuite/examples/associations.txt)


Run the command with the help option to see the usage and all the available options.

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

Using the tool to download all the sequences from the KO associated with the pathway with the ID map:

```bash
download_kos -o /path/to/output/folder/ -m map map00362
```

Using the tool to download all the sequences from the KO's listed in the file examples/kos.txt:

```bash
download_kos -o /path/to/output/folder/ -k /path/to/OrtScraper/examples/kos.txt
```

Using the tool to download all the sequences from the KO's associated with the listed reactions in the file examples/reactions.txt:

```bash
download_kos -o /path/to/output/folder/ -r /path/to/OrtScraper/examples/reactions.txt
```

Using the tool to download all the sequences from the KO's associated with the listed enzymes in the file examples/ecs.txt:

```bash
download_kos -o /path/to/output/folder/ -e /path/to/OrtScraper/examples/ecs.txt
```

Note: Running this commands may take some time and memory space.


Output
======
In the output folder you will find one FASTA file for each selected KO.
If you use one of the -e or -r options you will have another file, associations.txt, which indicates which kos where selected for download for each reaction/EC number.
In the same folder can be found a file info_db.csv that contains a table with information regarding the KO's that was selected to download, their name and also the EC numbers and Reactions IDs to which they are associated.

