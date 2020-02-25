# OrtAn

This tool uses [OrthoFinder](https://github.com/davidemms/OrthoFinder) output data to perform the annotation of the created Orthogroups based on a given database.

### Input:

- OrthoFinder working directory resulting from the run with the desired set of genomes.

- Database to base the annotation with.

### Dependencies:

Python 3

[DIAMOND](https://github.com/bbuchfink/diamond)

# Pipeline

Before start running OrtAn you need to prepare the input files:
- Run OrthoFinder with the input genomes;
- Prepare the database with the necessary format.

### OrthoFinder required files and directories:

```/Orthogroups/Orthogroups.txt```

```/WorkingDirectory/SpeciesIDs.txt```

```/WorkingDirectory/SequenceIDs.txt```

```/Orthogroup_Sequences/```


**Note:**  You don't need to run the entire workflow of OrthoFinder to get all the required files. You can use the option -og that stops after inferring the orthogroups.


### Database Format:

The database should be organized in a folder where each fasta file contain a set of (amino acid) sequences annotated for the same function.
You can create a database with this format using the [OrtScraper](https://github.com/MartaLopesGomes/OrtScraper) tool (Collect information from a given pathway or set of IDs from the KEGG database).


## Create Project

In this step, OrtAn receives the information of the input database and creates the working directory structure necessary for the following steps.

**Note:** you should indicate always the same output/working directory in all the steps.

Run ```create_project -h``` to see the usage of this command.

## Relaxed Search

This step can only be performed after you have run OrthoFinder with the desired genomes.

This command takes the OrthoFinder information and runs a first relaxed search with DIAMOND to identify the associations between the returned orthogroups and the different functions of the database.
On this first search, only 1 random sequence per 10 in each orthogroup is used to get the associations, reducing the search space on the next step.
Here you also can define the threshold of identity percent you want to use. Default: 80.

Run ```relaxed_search -h``` to see the usage of this command.

## Restrictive Search

In this step, OrtAn performs a restrictive search only between the orthogroups and groups of functions from the database that we find to be possibly related in the first step.

Run ```restrictive_search -h``` to see the usage of this command.


## Annotation

In this step, the tool uses the data created before to return the annotation of the sequences present in the orthogroups.
Here you can decide the thresholds of the parameters that will give you the annotation of the unknown sequences.

**% Identity** - the percentage of identical matches in the range of alignment. Default: 95.

**% Positive Matches** - percent of identical matches + positive matches in the range of alignment.Default: 99.

**% Query Coverage** - percent of the query sequence involved in the range of alignment.Default: 90.

**% Target Coverage** - percent of the target sequence (sequence in the database) involved in the range of alignment.Default: 90.

Run ```annotation -h``` to see the usage of this command.


## Create Database

With this command, you have the option to update your initial database adding the newly annotated sequences from the input genomes.
You can have the option to create a new database in a folder of your choice, or simply add the new sequences to the given database.

Run ```create_db -h``` to see the usage of this command.


## Instalation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install virtualenv.

```bash
pip install virtualenv
```

Create a virtual environment to use with this tool.

```bash
virtualenv orthoAnnotation
```

Activate the virtual environment.

```bash
source orthoAnnotation/bin/activate
```

Then move to the folder where the file setup.py from the OrtAn tool is located.

```bash
cd /path/to/OrtAn
```

Run the command:

```bash
python setup.py install
```

The tool should be ready to use.


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

From relaxed_search step, we get a text file (```/Results/Associations.txt```) containing the associations between the Orthogroups and the database functions.

From the annotation step we get 6 different text files:

```/Results/Annotation_Function_Protein.txt``` - Shows in the first column the functions and in the second the sequences annotated to that functions (one association per line).

```/Results/Annotation_Protein_Function.txt``` - Shows in the first column the sequences and the second the Functions for which the sequences were annotated for (one association per line).

```/Results/ConOG.txt``` - Consistent Orthogroups (Orthogroups where all the sequences were annotated to the same function). The function is also indicated.

```/Results/DivOG.txt``` - Divergent Orthogroups (Orthogroups where not all the sequences were annotated to the same function). This means that the orthogroup could have sequences that were not annotated to any function or sequences annotated to different functions. The functions are also indicated in the file.

```/Results/Orthogroups_Annotation.csv``` - This file shows how many sequences in each Orthogroup were annotated and to which function.

```/Results/Species_Annotation.csv``` - This file shows which functions are present in which species (1 - at least one sequence annotated to the function, 0 - no sequences annotated to the function).



