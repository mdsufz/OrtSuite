# OrtSuite - a flexible pipeline for annotation of ecosystem processes and prediction of putative microbial interactions

OrtSuite was developed with the goal to facilitate annotation of ecosystem processes and identify putative microbial interactions by automating the complete process from sequence retrieval, clustering of ortholog sequences, functional annotation, to putative microbial interactions prediction. 



![workflow](https://github.com/msdsufz/OrtSuite/blob/master/workflow_ortan_no_legend.png)
**OrtAn workflow** - Protein sequence from samples supplied by the user are clustered using OrthoFinder **(a)**. OrtScraper takes a text file containing a list of identifiers for each reaction in the pathway of interest supplied by the user to retrieve all protein sequences from KEGG **(b)**. Sequences mapped to reactions are stored in ORAdb **(c)**. Functional annotation **(d)** consists of a two-stage process (relaxed and restrictive search). Relaxed search **(e)** performs sequence alignments between 10% of randomly selected sequences from each generated cluster. Clusters whose representative sequences share a minimum 50% identity to sequences in reaction set(s) in ORAdb transition to the restrictive search **(f)**. Here, all sequences from the cluster is aligned to all sequences in the corresponding reaction set(s) to which they had a hit. Finally, the annotated sequences are used to identify putative microbial interactions based on their functional potential **(g and h)**. Additional constraints can be added to reduce the search space of microbial interactions **(rounded blue rectangles)**.


# Overview of OrtSuite

**OrtScraper:** Bulk download of protein sequences for populating a user-defined database

**OrthoFinder:** Clustering of orthologs

**OrtAn:** Functional annotation of clusters of orthologs


# System Requirements

Resources for Ortsuite will vary depending on the amount of data being processed. Nevertheless, we recommend a minimum of 4 cores and 8GB of RAM for small datasets (<50 identifiers). OrtSuite officially supports only Linux OS. 



Installation
============


# OrtScraper

OrtScraper is a python tool used to request information in bulk from KEGG (Kyoto Encyclopedia of Genes and Genomes) database.

**Requirements:**  Python 3

**Dependencies:**  bs4, grequests


We suggest the use of the package manager [pip](https://pip.pypa.io/en/stable/) to install the virtual environment (*virtualenv*).


```bash
pip install virtualenv
```

Create a virtual environment.

```bash
virtualenv venv_OrtScraper
```

Activate the virtual environment.

```bash
source venv_OrtScraper/bin/activate
```

Install dependencies

```bash
pip install grequests
pip install bs4
```

Next move to the folder where the file setup.py from the OrtScraper tool is located.

```bash
cd /path/to/OrtScraper
```

Run the command:

```bash
python setup.py install
```

Once the installation is finished the tool should be ready to use.



Usage
=====

## download_kos

Download all the sequences from each desired KO (KEGG Orthology) group to a FASTA file.

The input can be:

- KEGG pathway map ID

- List of KO IDs

- List of KEGG Reaction IDs

- List of EC (Enzyme commission) numbers

Note: The format of the input lists must be a txt file with only one ID per line.


Run the command with the help option to see the usage and all the available options.

```bash
download_kos -h
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

