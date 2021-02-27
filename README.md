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

#### data: Includes all the input data such as mutual_exclusivity files, MLA files, PPI files etc. "{method}_mutation_filtered_ep_data" folder contains the MEX files.

	A. binary_matrices_all_genes_ep_mutation_filtered: includes binary matrices for mutation data.

	B. {method}_mutation_filtered_ep_data: inlcudes pairwise mutual exclusivity p-values.

	C. MLA_ep_mutation_filtered_all_genes: Corresponding MLA.

	D: Census_allFri_Apr_26_12_49_57_2019.tsv: COSMIC file (Known driver genes)
	
	E: intact_nodupl_edge_file: intact network edge file 
	
	F: intact_nodupl_index_file: intact network index file 



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
