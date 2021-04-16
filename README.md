# A Network-centric Framework for the Evaluation of Mutual Exclusivity Tests on Cancer Drivers

This is the original repository for the Network-centric Framework for the Evaluation of Mutual Exclusivity Tests project codes. The project involves the evaluation of mutual exclusivity methods: Discover, Discover_strat, Fisher's Exact Test, MEGSA, MEMO, and WExT. The results from these methods include pairwise mutual exclusivity p-values. Based on them, we apply our network-centric epistatic evaluation.


## Getting Started
### Prerequisites
Pyhton: 3.8-3.9

### Getting Started

Using Github clone

```bash
git clone  https://github.com/abu-compbio/NetCentric 
cd NetCentric
pip install -r requirements.txt
```
You can also run this project locally by following these steps:
1. Download the repo
2. Unzip NetCentric-main
3. Open cmd/terminal and cd into the project
4. Execute python -m pip install -r requirements.txt

Some of the datasets are given in the .zip file format. In order to unzip, you should run script_unzip_data.py located under the NetCentric.

```bash
python script_unzip_data.py
```

## Input

1. The PPI network edge and index files. 

The file is located at data/intact_nodupl_index_file.txt 

```bash
0	""CHEBI
1	100147744
2	1B
3	1EFV
...
``` 
The file is located at data/intact_nodupl_edge_file.txt

```bash
7589	13441	0.99
9123	10446	0.98
4248	1740	0.98
3776	5279	0.98
...
``` 

2. Mutual Exclusivity Data

The mutation data includes pairwise mutual exclusivity p-values given for each method (discover, discover_strat, fishers, megsa, memo and wext).
The files with the name of mutations_all_genes include all genes and intact_filtered include only ones in intact network. 

The file is located at data/{method}_mutation_filtered_ep_data/{cancerType}_{method}_result_mutations_all_genes_{threshold}.txt
```bash
	gene1	gene2	pvalue	qvalue
0	A2M	A2ML1	0.6584654889330113	0.9926886078786901
1	A2M	ABCA1	0.5332913581418495	0.9926886078786901
2	A2M	ABCA10	0.8971732886956303	0.9926886078786901
...
``` 
The file is located at data/{method}_mutation_filtered_ep_data/{cancerType}_{method}_pairs_intact_filtered_subset{threshold}.txt

```bash
	gene1	gene2	pvalue	oddsratio
0	TCF7L2	CTNNB1	0.9015805073650888	1.6786858974358974
1	SMAD4	SMAD3	0.839665475908354	1.4567901234567902
2	EP300	TP53	0.0742406168447767	0.5221052631578947
...
``` 

Binary genes matrix: The file contains binary martices for mutation data. 

The file is located at data/binary_matrices_all_genes_ep_mutation_filtered

```bash
	patients	A1BG	A1CF	A2M ...
TCGA-3L-AA1B-01A	0	0	0
TCGA-4N-A93T-01A	0	0	0
TCGA-4T-AA8H-01A	0	0	0
...
``` 

MLA: The file contains the corresponding MLA.

The file is located at data/MLA_ep_mutation_filtered_all_genes

```bash
A1BG	4.261253658699028
A1CF	5.095042391780406
A2M	5.539871662596874
...
``` 

3. COSMIC File (Driver genes)

The file is located at data/Census_allFri_Apr_26_12_49_57_2019.tsv

4. TSN

The file contains gtex edges for the corresponding tissue type with a given threshold (0.0 and 0.5)

The file is located at data/gtex_tsn_fractions_intact_filtered_applied_threshold

```bash
MDM2	TP53	1.0
PAK1	RAC1	1.0
FADD	CASP8	0.9987163029525032
...
``` 

## Runs

The codes regarding various analyses given in the main article.

### **ME Evaluations Based on Defined Metrics** 

The main source code for the evaluation ME Tests. As output, you get tables with all analysis results in NetCentric/results_main/evaluation_results
To generate the algorithm for the given input, the following script should be run

```bash
cd src
python evaluations_on_metrics.py
``` 

### **ME Evaluations Based on Corrections via MLA**

Scatterplots of percentage significance of mutual exclusivity runs vs mutation load association (MLA)
As output, you get results in NetCentric/results_main/evaluation_results/percent_sig_figures

```bash
cd src
python evaluations_via_mla.py
```

Scatterplots of percentage significance of mutual exclusivity runs vs mutation load association (MLA) when only CGC genes that have > 1 neighbors are included.
As output, you get results in NetCentric/results_main/evaluation_results/percent_sig_figures/perc_sig_figures_for_multiple_neighbors

```bash
cd src
python evaluations_via_mla_neighbors.py
```

### **ME Evaluations Based on Corrections via TSN**

Network-centric epistasis evaluation framework run on the tissue-specific network (TSN). As output, you get results in NetCentric/results_main/evaluation_results/intact

```bash
cd src
python evaluations_via_tsn.py
```

ROC curves for comparing the mutual exclusivities of tissue-specific and non-tissue specific CGC-CGC gene pairs and non-CGC-non-CGC gene pairs on.
As output, you get results in NetCentric/results_main/figure_tsn_AUROC

```bash
cd src
python me_on_tsn_ntsn_roc_curve.py

```

## Output

The outputs;

* The evaluation results for the IntAct network and TSN network,
* Scatterplots of percentage significance of mutual exclusivity runs vs mutation load association (MLA),
* ROC curves for comparing the mutual exclusivities of tissue-specific and non-tissue-specific CGC-CGC gene pairs and non-CGC-non-CGC gene pairs

will be stored in the results_main file. 

