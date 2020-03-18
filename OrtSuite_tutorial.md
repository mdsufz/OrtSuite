Guide to predict putative microbial interactions with OrtSuite
====

Briefly, OrtSuite consists of **1:** Generating the Ortholog Reaction-Association database (*ORAdb*), **2:** Clustering of Orthologs, **3:** Functional annotation of clusters of orthologs, and **4:** Identification of putative microbial interactions.


1: Generating the Ortholog Reaction-Association database
====

Download all sequences associated with the list of enzyme commission numbers:
```bash

First make a directory: 
mkdir test_Ortsuite

download_kos -o ~/test_OrtSuite/ -e /OrtSuite/examples/ecs.txt

```
