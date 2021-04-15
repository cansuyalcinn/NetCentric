# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 22:59:48 2021

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm,trange
import operator
from scipy.stats import mannwhitneyu
from sklearn.metrics import auc
import random
random.seed(1234)
import os
import string


dict_tissue={'BLCA': 'Bladder',
            'BRCA': 'Breast',
            'COADREAD': 'Colon',
            'LUAD': 'Lung',
            'LUSC': 'Lung',
            'SKCM': 'Skin',
            'STAD': 'Stomach',
            'UCEC': 'Uterus'
            }


dict_methods = {'discover': 'DISCOVER',
               'discover_strat': 'DISCOVER Strat',
               'fishers': 'Fisher\'s Exact Test',
               'megsa': 'MEGSA',
               'memo': 'MEMO',
               'wext': 'WExT'}
#mutation threshold: t, cancer type: c, tsn threshold: threshold, perc=percentage
c = 'COADREAD'
methods = ['discover','discover_strat','fishers','megsa', 'memo','wext']
t=20
tissue = dict_tissue[c]
threshold = 0.0
perc = 0.25

#inputs
inpath_mla = '../data/MLA_ep_mutation_filtered_all_genes/'
cosmic_infile = '../data/Census_allFri_Apr_26_12_49_57_2019.tsv'
intact_edge_file = '../data/intact_nodupl_edge_file.txt'
intact_index_file = '../data/intact_nodupl_index_file.txt'
infile_tsn = '../data/gtex_tsn_fractions_intact_filtered_applied_threshold/edges_gtex_intact_filtered_{}_{}.txt'.format(tissue,threshold)

with open(infile_tsn) as f:
    tsn_edges = {(line.split()[0].upper(),line.split()[1].upper()):float(line.split()[2]) for line in f.readlines()}
    
with open(cosmic_infile,'r') as f:
    cosmic_genes = [line.split()[0].upper() for line in f.readlines()[1:]]

with open(intact_index_file, 'r') as f:
    indices = {line.split()[0]:line.split()[1] for line in f.readlines()}

with open(intact_edge_file, 'r') as f:
    edges = [(indices[line.split()[0]].upper(),indices[line.split()[1]].upper()) for line in f.readlines()]

intact_genes_list = list(indices.values())

count = 0
with open(infile_tsn) as f:
    lines_tsn=f.readlines()
    dict_tsn_conf = {(line.split()[0],line.split()[1]):float(line.split()[2]) for line in lines_tsn}
    tsn_top_edges = [(line.split()[0].upper(),line.split()[1].upper()) for line in lines_tsn if float(line.split()[2])>=(1-perc)]
    tsn_bottom_edges = [(line.split()[0].upper(),line.split()[1].upper()) for line in lines_tsn if float(line.split()[2])<=perc]
    
    
print(len(tsn_edges), len(dict_tsn_conf),len(tsn_top_edges), len(tsn_bottom_edges))

inpath = '../data/binary_matrices_all_genes_ep_mutation_filtered/'

df = pd.read_csv(inpath + c + '_TML_binary_sm.txt', sep= '\t', index_col=0, header=0, skipinitialspace=True)

df.drop('y',1,inplace=True)

df.drop([col for col, val in df.sum().iteritems() if val <= t], axis=1, inplace=True)

dict_inpath = {'discover': '../data/discover_mutation_filtered_ep_data/',
               'discover_strat': '../data/discover_mutation_filtered_ep_data',
               'fishers': '../data/fishers_mutation_filtered_ep_data/',
               'megsa': '../data/megsa_mutation_filtered_ep_data/',
               'memo': '../data/memo_mutation_filtered_ep_data/',
               'wext': '../data/wext_mutation_filtered_ep_data/'
              }

inpath_mla = '../data/MLA_ep_mutation_filtered_all_genes/'
cosmic_infile = '../data/Census_allFri_Apr_26_12_49_57_2019.tsv'

dict_infile = {}
dict_infile_intact = {}
for m in methods:
    if m == 'discover':
            suffix = '/{}_{}_result_mutations_all_genes_q1.0_normal_{}.txt'.format(c,m,t)
            suffix_intact = '/{}_pairs_q1.0_normal_intact_filtered_subset{}.txt'.format(c,t)
    elif m == 'discover_strat':
            suffix = '/{}_{}_result_mutations_all_genes_q1.0_stratified_{}.txt'.format(c,'discover',t)
            suffix_intact = '/{}_pairs_q1.0_stratified_intact_filtered_subset{}.txt'.format(c,t)    
    else:
            suffix = '/{}_{}_result_mutations_all_genes_{}.txt'.format(c,m,t)
            suffix_intact = '/{}_{}_pairs_intact_filtered_subset{}.txt'.format(c,m,t)

    dict_infile[m] = dict_inpath[m] + suffix
    dict_infile_intact[m] = dict_inpath[m] + suffix_intact 

def get_pvalues(filename, reference_genes=cosmic_genes,pvalue_threshold=0.05,pvalue_position=3, ref_edges=tsn_edges):
    """ge mutex pvalues from file"""
          
    with open(filename, 'r') as f:
        lines = f.readlines()[1:]
        
        list_pvals = []
        for line in tqdm(lines):
            val = float(line.split()[pvalue_position])
            if val!=0:
                list_pvals.append(val)
        min_pval = min(list_pvals)
        del list_pvals
    print('min pval:',min_pval)

    dict_temp = {g:{} for g in tsn_genes}

    for line in tqdm(lines):
        line = line.strip().split('\t')

        g1 = line[1]
        g2 = line[2]

        p_temp = float(line[pvalue_position])

        if p_temp == 0:
            if m=='wext':
                p=0.0
            else:
                p = -np.log(min_pval)
        else:
            p = -np.log(p_temp)

        if g1 in dict_temp and g2 not in dict_temp[g1]:
            dict_temp[g1][g2] = p
        if g2 in dict_temp and g1 not in dict_temp[g2]:
            dict_temp[g2][g1] = p

    return dict_temp
    

tsn_genes = set()
for g1,g2 in tsn_edges:
    tsn_genes.update([g1,g2])
tsn_genes.intersection_update(set(df.columns.tolist()))

dict_mex = {}
for m in tqdm(methods):
    print(m)
    filename = dict_infile[m]
    dict_mex[m]=get_pvalues(filename)


genes = df.columns.to_list()
dict_tsn_conf_cosmic_nb = {k:v for k,v in dict_tsn_conf.items() if k[0] in cosmic_genes and k[1] in cosmic_genes and k[0] in genes and k[1] in genes}

top_n=1.0
bottom_n = 0.5
dict_tsn_conf_cosmic_nb_top = {k:v for k,v in dict_tsn_conf_cosmic_nb.items() if v>=top_n}
dict_tsn_conf_cosmic_nb_bottom = {k:v for k,v in dict_tsn_conf_cosmic_nb.items() if v<=bottom_n}
len(dict_tsn_conf_cosmic_nb_top),len(dict_tsn_conf_cosmic_nb_bottom)

outpath = '../results_main/figure_tsn_AUROC/'

if not os.path.exists(outpath):
    os.makedirs(outpath)
    
fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(30,16))
axes = np.array(axes)
axiter = axes.flat

for i, (ax,m) in enumerate(zip(axiter,methods)):

    # COSMIC-COSMIC
    list_vals = []

    for k,v in dict_tsn_conf.items():
        g1,g2=k
        if g1 in cosmic_genes and g2 in cosmic_genes and g1 in genes and g2 in genes: 
            if m!='wext':
                if v>=top_n:
                    list_vals.append([g1,g2,dict_mex[m][g1][g2],1])
                elif v<=bottom_n:
                    list_vals.append([g1,g2,dict_mex[m][g1][g2],0])
            else:
                if v>=top_n:
                    if g2 in dict_mex[m][g1]:
                        list_vals.append([g1,g2,dict_mex[m][g1][g2],1])
                    else:
                        continue
                elif v<=bottom_n:
                    if g2 in dict_mex[m][g1]:
                        list_vals.append([g1,g2,dict_mex[m][g1][g2],0])
                    else:
                        continue
            
            
    df_rank = pd.DataFrame(list_vals, columns=['g1', 'g2', 'mex','class'])
    df_rank.sort_values('mex', inplace=True,ascending=False)

    '''
    TPR = Sensitivity = sum_true/len(dict_tsn_conf_cosmic_nb_top),
    FPR = 1-Specificity = Sum_false/len(dict_tsn_conf_cosmic_nb_bottom).
    '''

    list_sum_true =[]
    list_sum_false = []

    count_true=0
    count_false=0
    for idx, row in df_rank.iterrows():
        cls = row['class']
        if cls==0:
            count_false+=1
        elif cls == 1:
            count_true+=1
        list_sum_true.append(count_true)
        list_sum_false.append(count_false)

    tpfn = float(len(dict_tsn_conf_cosmic_nb_top))
    fptn = float(len(dict_tsn_conf_cosmic_nb_bottom))
    df_rank['sum_true'] = list_sum_true
    df_rank['sum_false'] = list_sum_false
    df_rank['tpr'] = [float(v)/tpfn for v in list_sum_true]
    df_rank['fpr'] = [float(v)/fptn for v in list_sum_false]
      
    ax.plot(df_rank['fpr'].tolist(),df_rank['tpr'].tolist(),c='C0', linewidth=4.5)

    ax.text(-0.1, 1.1, string.ascii_uppercase[i], transform=ax.transAxes, 
            size=25, weight='bold')
    
    #Noncosmic-Noncosmic
    list_auc=[]

    list_vals2_temp_top = []
    list_vals2_temp_bottom = []
    list_vals2 = []
    
    for k,v in dict_tsn_conf.items():
        g1,g2=k
        if g1 not in cosmic_genes and g2 not in cosmic_genes and g1 in genes and g2 in genes: 
            if m!='wext':
                if v>=top_n:
                    list_vals2_temp_top.append([g1,g2,dict_mex[m][g1][g2],1])
                elif v<=bottom_n:
                    list_vals2_temp_bottom.append([g1,g2,dict_mex[m][g1][g2],0])
            else:
                if v>=top_n:
                    if g2 in dict_mex[m][g1]:
                        list_vals2_temp_top.append([g1,g2,dict_mex[m][g1][g2],1])
                    else:
                        list_vals2_temp_top.append([g1,g2,0,1])    
                elif v<=bottom_n:
                    if g2 in dict_mex[m][g1]:
                        list_vals2_temp_bottom.append([g1,g2,dict_mex[m][g1][g2],0])
                    else:
                        continue
    list_tpr = []
    list_fpr = []
    for i in trange(100):

        list_vals2 = random.sample(list_vals2_temp_top,len(dict_tsn_conf_cosmic_nb_top))+random.sample(list_vals2_temp_bottom,len(dict_tsn_conf_cosmic_nb_bottom))

        df_rank_nc = pd.DataFrame(list_vals2, columns=['g1', 'g2', 'mex','class'])
        df_rank_nc.sort_values('mex', inplace=True,ascending=False)

        list_sum_true_nc =[]
        list_sum_false_nc = []

        count_true_nc=0
        count_false_nc=0
        for idx, row in df_rank_nc.iterrows():
            cls = row['class']
            if cls==0:
                count_false_nc+=1
            elif cls == 1:
                count_true_nc+=1
            list_sum_true_nc.append(count_true_nc)
            list_sum_false_nc.append(count_false_nc)

        tpfn_nc = float(count_true_nc)
        fptn_nc = float(count_false_nc)
        df_rank_nc['sum_true'] = list_sum_true_nc
        df_rank_nc['sum_false'] = list_sum_false_nc
        df_rank_nc['tpr'] = [float(v)/tpfn_nc for v in list_sum_true_nc]
        df_rank_nc['fpr'] = [float(v)/fptn_nc for v in list_sum_false_nc]

        auc_score = auc(df_rank_nc['fpr'],df_rank_nc['tpr'])
        list_auc.append(auc_score)
        list_fpr.append(df_rank_nc['fpr'].tolist())
        list_tpr.append(df_rank_nc['tpr'].tolist())

    med_auc_idx = np.argsort(list_auc)[len(list_auc)//2]
    med_auc = list_auc[med_auc_idx]

    ax.plot(list_fpr[med_auc_idx],list_tpr[med_auc_idx], c='C1', linewidth=4.5)
    
    ax.set_xlabel('False Positive Rate',fontsize=20) 
    ax.set_ylabel('True Positive Rate',fontsize=20)
    ax.set_title(dict_methods[m],fontsize=25)
    for tick in ax.xaxis.get_major_ticks(): 
        tick.label.set_fontsize(20)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)

    legend_ = []

    legend_.append('{} ({:.2f})'.format('CGC-CGC pairs', auc(df_rank['fpr'],df_rank['tpr'])))
    legend_.append('{} ({:.2f})'.format('nonCGC-nonCGC pairs', med_auc))


    art = []
    legend = ax.legend(legend_, loc=8,fancybox=True, fontsize=18, framealpha=0,
                        edgecolor = 'b', ncol= 2, bbox_to_anchor=(0.5,-0.2))
    art.append(legend)
    frame = legend.get_frame()

fig.tight_layout() 
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Calibri' 
plt.savefig(outpath+'{}_t{}_tsn_auroc.pdf'.format(c,t),format='pdf', bbox_inches='tight')

