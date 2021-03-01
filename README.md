# A Network-centric Framework for the Evaluation of Mutual Exclusivity Tests on Cancer Drivers

This is the original repository for the Network-centric Framework for the Evaluation of Mutual Exclusivity Tests project codes. The project involves the evaluation of mutual exclusivity methods: Discover, Discover_strat, Fisher's Exact Test, MEGSA, MEMO, and WExT. The results from these methods include pairwise mutual exclusivity p-values. Based on them, we apply our network-centric epistatic evaluation.


## Installation

Installing the dependencies

```bash
conda install --file requirements.txt
```
Installing the conda env

```bash
conda env create -f netcent2.yml
``` 

## Input

#### data
1. The PPI network edge and index files. 

The file is located at data/intact_nodupl_edge_file.txt

```bash
0	""CHEBI
1	100147744
2	1B
3	1EFV
...
``` 
The file is located at data/intact_nodupl_index_file.txt 

```bash
7589	13441	0.99
9123	10446	0.98
4248	1740	0.98
3776	5279	0.98
...
``` 

2. Mutation data

The mutation data includes pairwise mutual exclusivity p-values given for each method (discover, discover_strat, fishers, megsa, memo and wext).
The files with the name of mutations_all_genes include all genes and intact_filtered include only ones in intact network. 

The file is located at data/{method}_mutation_filtered_ep_data/{cancer type}_{method}_result_mutations_all_genes_{threshold}

```bash
	gene1	gene2	pvalue	qvalue
0	A2M	A2ML1	0.6584654889330113	0.9926886078786901
1	A2M	ABCA1	0.5332913581418495	0.9926886078786901
2	A2M	ABCA10	0.8971732886956303	0.9926886078786901
...
``` 
The file is located at data/{method}_mutation_filtered_ep_data/{cancer type}_{method}_pairs_intact_filtered_subset{threshold}

```bash
	gene1	gene2	pvalue	oddsratio
0	TCF7L2	CTNNB1	0.9015805073650888	1.6786858974358974
1	SMAD4	SMAD3	0.839665475908354	1.4567901234567902
2	EP300	TP53	0.0742406168447767	0.5221052631578947
...
``` 

Binary genes matrix:
The file contains binary martices for mutation data. 

The file is located at data/binary_matrices_all_genes_ep_mutation_filtered

```bash
	patients	A1BG	A1CF	A2M ...
TCGA-3L-AA1B-01A	0	0	0
TCGA-4N-A93T-01A	0	0	0
TCGA-4T-AA8H-01A	0	0	0
...
``` 

2. COSMIC File (Driver genes)

The file is located at data/Census_allFri_Apr_26_12_49_57_2019.tsv





#### data: Includes all the input data such as mutual_exclusivity files, MLA files, PPI files etc. "{method}_mutation_filtered_ep_data" folder contains the MEX files.

	A. binary_matrices_all_genes_ep_mutation_filtered: includes binary matrices for mutation data.

	C. MLA_ep_mutation_filtered_all_genes: Corresponding MLA.
	
	G: gtex_tsn_fractions_intact_filtered_applied_threshold: gtex edges threshold 0.0 and 0.5 were applied



## Runs

The codes regarding various analyses given in the main article.

**ME Evaluations Based on Defined Metrics** 

The main source code for the evaluationon ME Tests. As output, you get tables with all analysis results.

* netcentric_me_evaluations_based_on_defined_metrics.py

**ME Evaluations Based on Corrections via MLA**

Scatterplots of percentage significance of mutual exclusivity runs vs mutation load
association (MLA)

* me_evaluations_based_on_corrections_via_MLA.py

Scatterplots of percentage significance of mutual exclusivity runs vs mutation load
association (MLA) when only CGC genes that have > 1 neighbors are included

* me_evaluations_based_on_corrections_via_MLA_neighbors.py

**ME Evaluations Based on Corrections via TSN**

Network-centric epistasis evaluation framework run on the tissue-specific network (TSN)

* me_evaluations_based_on_corrections_via_TSN.py

ROC curves for comparing the mutual exclusivities of tissue-specific and non-tissue specific CGC-CGC gene pairs and non-CGC-non-CGC gene pairs on

* me_on_tsn_and_ntsn_roc_curve.py
