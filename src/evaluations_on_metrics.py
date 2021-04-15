# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 21:15:02 2021

@author: Rafsan Ahmed & Cansu Yalcin
"""
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
import gc
from tqdm.notebook import tqdm
from scipy.stats import pearsonr,spearmanr
from random import sample,choice
import random 
random.seed(1234)

pd.set_option('display.max_columns', None)

#mutation threshold: t, cancer type: c
t = 20
c = 'COADREAD' 
methods = ['discover','discover_strat','fishers','megsa','memo','wext']

#inputs
cosmic_infile = '../data/Census_allFri_Apr_26_12_49_57_2019.tsv'
save_path = '../results_main/evaluation_results/intact'
intact_edge_file = '../data/intact_nodupl_edge_file.txt'
intact_index_file = '../data/intact_nodupl_index_file.txt'

dict_infile = {}
dict_infile_intact = {}
for m in methods:
    if m == 'discover':
            suffix = '{}_mutation_filtered_ep_data/{}_{}_result_mutations_all_genes_q1.0_normal_{}.txt'.format(m,c,m,t)
            suffix_intact = '{}_mutation_filtered_ep_data/{}_pairs_q1.0_normal_intact_filtered_subset{}.txt'.format(m,c,t)
    elif m == 'discover_strat':
            suffix = '{}_mutation_filtered_ep_data/{}_{}_result_mutations_all_genes_q1.0_stratified_{}.txt'.format('discover',c,'discover',t)
            suffix_intact = '{}_mutation_filtered_ep_data/{}_pairs_q1.0_stratified_intact_filtered_subset{}.txt'.format('discover',c,t)    
    else:
            suffix = '{}_mutation_filtered_ep_data/{}_{}_result_mutations_all_genes_{}.txt'.format(m,c,m,t)
            suffix_intact = '{}_mutation_filtered_ep_data/{}_{}_pairs_intact_filtered_subset{}.txt'.format(m,c,m,t)

    dict_infile[m] = '../data/' + suffix
    dict_infile_intact[m] = '../data/' + suffix_intact


with open(intact_index_file, 'r') as f:
    indices = {line.split()[0]:line.split()[1] for line in f.readlines()}

with open(intact_edge_file, 'r') as f:
    edges = [(indices[line.split()[0]].upper(),indices[line.split()[1]].upper()) for line in f.readlines()]
    
with open(cosmic_infile,'r') as f:
    cosmic_genes = [line.split()[0].upper() for line in f.readlines()[1:]]
    
def chunks(list_of_genes,n=1000):
    """Seperate total genes into chunks for memory management"""
    for i in range(0,len(list_of_genes),n):
        yield list_of_genes[i:i+n]
        
def get_genes(filename):
    """get all genes in cohort, return a list"""
    with open(filename, 'r') as f:
        genes = set()
        for line in tqdm(f.readlines()[1:],desc='Counting total Genes'):
            genes.update(line.strip().split('\t')[1:3])

    return list(genes)   

def get_neighbors(genes, ref_edges):
    """create a dictionary of neighbors of genes. for each gene in the dict, all its neighbors
    will be present in the subdictionary.
    ref_edges = g1,g2 must be a part of the reference edges (ppi network edges)"""
    dict_neighbors = {}
    for g1,g2 in ref_edges:
        if g1 in genes and g2 in genes:
            if g1 not in dict_neighbors:
                dict_neighbors[g1] = set()
            if g2 not in dict_neighbors:
                dict_neighbors[g2] = set()

            dict_neighbors[g1].update([g2])
            dict_neighbors[g2].update([g1])

    return dict_neighbors

def get_cg_cg_genes(cohort_specific_genes, dict_neighbor,ref_genes=cosmic_genes):
    """get (cosmic gene --- cosmic gene) pairs"""
    set_cg_cg = set()
    for g in set.intersection(set(cohort_specific_genes),set(dict_neighbor), set(ref_genes)):
        if len(set.intersection(set(dict_neighbor[g]), set(ref_genes)))>0:
               set_cg_cg.update([g])
                
    return set_cg_cg

    
def count_sig_cosmic_pairs(filename, reference_genes=cosmic_genes,sig_threshold = 0.05):
    '''get the count of significant cosmic --- cosmic pairs'''
    count=0
    with open(filename, 'r') as f:
        for line in tqdm(f.readlines()[1:]):
            line=line.split()
            g1,g2,p = line[1],line[2],float(line[3])
            
            if (g1 in reference_genes or g2 in reference_genes) and p<sig_threshold:
                count+=1
    del line
    gc.collect()
    
    return count
    
def get_sig_logpval_counts_cgcg_minus_cgnnb_single(d, neighbor_set,reference_genes=cosmic_genes,randiter=100,sig_threshold=-np.log(0.05),zero_threshold=0):
    '''
    First evaluation: for each CGC gene g, check all its neighbors. 
    If more CGC-CGC pairs are present, get the necessary statistics by comparing CGC-CGC pairs vs CGC-non neighbor pairs. 
    d: dictionary for single CGC g containing all its neighbors,
    neighbor_set: all neighbors of g,
    randiter: inner iteration.
    
    '''
    d_cg_cg = {k:v for k,v in d.items() if k in reference_genes and k in neighbor_set}
    d_cg_nnb = {k:v for k,v in d.items() if k in reference_genes and k not in neighbor_set}
    
    if len(d_cg_cg)==0 or len(d_cg_cg)>len(d_cg_nnb):
        return (np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)

    else:
        count_sigLHS_nonsigRHS = []
        count_sigLHS_sigRHS_LHS = [] 
        count_sigLHS_sigRHS_RHS = [] 
        count_nonsigLHS_sigRHS = [] 
        count_nonsigLHS_nonsigRHS_LHS = [] 
        count_nonsigLHS_nonsigRHS_RHS = [] 
        sig_sum_RHS = []
        sum_RHS = []
        
        sig_sum_LHS = np.sum([v for v in d_cg_cg.values() if v>sig_threshold])
        sum_LHS = np.sum([v for v in d_cg_cg.values()])
        
        for i in range(randiter):
            count_sigLHS_nonsigRHS_temp = 0 
            count_sigLHS_sigRHS_LHS_temp = 0 
            count_sigLHS_sigRHS_RHS_temp = 0 
            count_nonsigLHS_sigRHS_temp = 0 
            count_nonsigLHS_nonsigRHS_LHS_temp = 0 
            count_nonsigLHS_nonsigRHS_RHS_temp = 0 
            sig_sum_RHS_temp = 0
            sum_RHS_temp = 0
            
            for cg in d_cg_cg:
                
                pval_cgcg = d_cg_cg[cg]
                
                rand_cg_nnb =choice(list(d_cg_nnb))
                pval_rand_cg_nnb = d_cg_nnb[rand_cg_nnb]
                sum_RHS_temp+=pval_rand_cg_nnb
            
                if pval_cgcg > sig_threshold:
                    
                    if pval_rand_cg_nnb < sig_threshold:
                        count_sigLHS_nonsigRHS_temp+=1
                    else:
                        sig_sum_RHS_temp+=pval_rand_cg_nnb
                        
                        if pval_cgcg > pval_rand_cg_nnb:
                            count_sigLHS_sigRHS_LHS_temp+=1
                        else:
                            count_sigLHS_sigRHS_RHS_temp+=1
                            
                else: 
                    if pval_rand_cg_nnb > sig_threshold:
                        sig_sum_RHS_temp+=pval_rand_cg_nnb
                        count_nonsigLHS_sigRHS_temp+=1
                    else:
                        
                        if pval_cgcg > pval_rand_cg_nnb:
                            count_nonsigLHS_nonsigRHS_LHS_temp+=1
                        else:
                            count_nonsigLHS_nonsigRHS_RHS_temp+=1         
                        
            count_sigLHS_nonsigRHS.append(count_sigLHS_nonsigRHS_temp) 
            count_sigLHS_sigRHS_LHS.append(count_sigLHS_sigRHS_LHS_temp) 
            count_sigLHS_sigRHS_RHS.append(count_sigLHS_sigRHS_RHS_temp) 
            count_nonsigLHS_sigRHS.append(count_nonsigLHS_sigRHS_temp) 
            count_nonsigLHS_nonsigRHS_LHS.append(count_nonsigLHS_nonsigRHS_LHS_temp) 
            count_nonsigLHS_nonsigRHS_RHS.append(count_nonsigLHS_nonsigRHS_RHS_temp)
            sig_sum_RHS.append(sig_sum_RHS_temp)
            sum_RHS.append(sum_RHS_temp)
        
        return len(d_cg_cg),np.median(count_sigLHS_nonsigRHS), np.median(count_sigLHS_sigRHS_LHS), np.median(count_sigLHS_sigRHS_RHS),\
        np.median(count_nonsigLHS_nonsigRHS_LHS),np.median(count_nonsigLHS_nonsigRHS_RHS), np.median(count_nonsigLHS_sigRHS),\
        sum_LHS,sig_sum_LHS, np.median(sum_RHS),np.median(sig_sum_RHS),\
        np.median(count_sigLHS_nonsigRHS)/float(len(neighbor_set)),\
        np.median(count_sigLHS_sigRHS_LHS)/float(len(neighbor_set)),\
        np.median(count_sigLHS_sigRHS_RHS)/float(len(neighbor_set)),\
        np.median(count_nonsigLHS_nonsigRHS_LHS)/float(len(neighbor_set)),\
        np.median(count_nonsigLHS_nonsigRHS_RHS)/float(len(neighbor_set)),\
        np.median(count_nonsigLHS_sigRHS)/float(len(neighbor_set))


def get_sig_logpval_counts_cgcg_minus_cgncgnb_single(d, neighbor_set,reference_genes=cosmic_genes,randiter=100,sig_threshold=-np.log(0.05),zero_threshold=0):
    '''
    Second evaluation: for each CGC gene g, check all its neighbors. 
    If more CGC-CGC pairs are present, get the necessary statistics by comparing CGC-CGC pairs vs CGC-non CGC neighbor pairs. 
    d: dictionary for single CGC g containing all its neighbors,
    neighbor_set: all neighbors of g,
    randiter: inner iteration,
    sig_threshold: significance threshold for MEX values, commonly 0.05.
    '''
    d_cg_cg = {k:v for k,v in d.items() if k in reference_genes and k in neighbor_set}
    d_cg_ncgnb = {k:v for k,v in d.items() if k not in reference_genes and k in neighbor_set}
    
    if len(d_cg_cg)==0 or len(d_cg_cg)>len(d_cg_ncgnb):
        return (np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan)

    else:
        count_sigLHS_nonsigRHS = [] 
        count_sigLHS_sigRHS_LHS = [] 
        count_sigLHS_sigRHS_RHS = [] 
        count_nonsigLHS_sigRHS = [] 
        count_nonsigLHS_nonsigRHS_LHS = [] 
        count_nonsigLHS_nonsigRHS_RHS = [] 
        sig_sum_RHS = []
        sum_RHS = []
        
        sig_sum_LHS = np.sum([v for v in d_cg_cg.values() if v>sig_threshold])
        sum_LHS = np.sum([v for v in d_cg_cg.values()])
        

        for i in range(randiter):
            count_sigLHS_nonsigRHS_temp = 0 
            count_sigLHS_sigRHS_LHS_temp = 0 
            count_sigLHS_sigRHS_RHS_temp = 0 
            count_nonsigLHS_sigRHS_temp = 0 
            count_nonsigLHS_nonsigRHS_LHS_temp = 0 
            count_nonsigLHS_nonsigRHS_RHS_temp = 0 
            sig_sum_RHS_temp = 0
            sum_RHS_temp = 0
            
            for cg in d_cg_cg:
                
                pval_cgcg = d_cg_cg[cg]
                
                rand_cg_ncgnb =choice(list(d_cg_ncgnb))
                pval_rand_cg_ncgnb = d_cg_ncgnb[rand_cg_ncgnb]
                sum_RHS_temp+=pval_rand_cg_ncgnb
                

                if pval_cgcg > sig_threshold:
                    
                    if pval_rand_cg_ncgnb < sig_threshold:
                        count_sigLHS_nonsigRHS_temp+=1
                    else:
                        sig_sum_RHS_temp+=pval_rand_cg_ncgnb
                        
                        if pval_cgcg > pval_rand_cg_ncgnb:
                            count_sigLHS_sigRHS_LHS_temp+=1
                        else:
                            count_sigLHS_sigRHS_RHS_temp+=1
                            
                else: 
                    if pval_rand_cg_ncgnb > sig_threshold:
                        sig_sum_RHS_temp+=pval_rand_cg_ncgnb
                        count_nonsigLHS_sigRHS_temp+=1
                    else:
                        
                        if pval_cgcg > pval_rand_cg_ncgnb:
                            count_nonsigLHS_nonsigRHS_LHS_temp+=1
                        else:
                            count_nonsigLHS_nonsigRHS_RHS_temp+=1         
                        
            count_sigLHS_nonsigRHS.append(count_sigLHS_nonsigRHS_temp) 
            count_sigLHS_sigRHS_LHS.append(count_sigLHS_sigRHS_LHS_temp) 
            count_sigLHS_sigRHS_RHS.append(count_sigLHS_sigRHS_RHS_temp) 
            count_nonsigLHS_sigRHS.append(count_nonsigLHS_sigRHS_temp) 
            count_nonsigLHS_nonsigRHS_LHS.append(count_nonsigLHS_nonsigRHS_LHS_temp)
            count_nonsigLHS_nonsigRHS_RHS.append(count_nonsigLHS_nonsigRHS_RHS_temp)
            sig_sum_RHS.append(sig_sum_RHS_temp)
            sum_RHS.append(sum_RHS_temp)
        
        
        return len(d_cg_cg),np.median(count_sigLHS_nonsigRHS), np.median(count_sigLHS_sigRHS_LHS), np.median(count_sigLHS_sigRHS_RHS),\
        np.median(count_nonsigLHS_nonsigRHS_LHS),np.median(count_nonsigLHS_nonsigRHS_RHS), np.median(count_nonsigLHS_sigRHS),\
        sum_LHS,sig_sum_LHS, np.median(sum_RHS),np.median(sig_sum_RHS),\
        np.median(count_sigLHS_nonsigRHS)/float(len(neighbor_set)),\
        np.median(count_sigLHS_sigRHS_LHS)/float(len(neighbor_set)),\
        np.median(count_sigLHS_sigRHS_RHS)/float(len(neighbor_set)),\
        np.median(count_nonsigLHS_nonsigRHS_LHS)/float(len(neighbor_set)),\
        np.median(count_nonsigLHS_nonsigRHS_RHS)/float(len(neighbor_set)),\
        np.median(count_nonsigLHS_sigRHS)/float(len(neighbor_set))
    
def transaction(df):
    """
    Changing the format of the table in accordance with the article.
    df : Dataframe containing the result table.
    df1_new : The final dataframe after format conversion is done.
    """
    columns = df.columns.tolist()
    columns = columns[:1]+columns[2:8]+columns[12:14]+columns[-4:]
    new_cols = ['Method', 'Stat1','Stat2','Stat3','Stat4','Stat5','Stat6', 'Stat7', 'Stat8', 'Precision', 'Sensitivity', 'Specificity', 'F1 Score']
    dict_methods = {'discover': 'DISCOVER',
                   'discover_strat': 'DISCOVER Strat',
                   'fishers': 'Fisher\'s Exact Test',
                   'megsa': 'MEGSA',
                   'memo': 'MEMO',
                   'wext': 'WExT',}
    decimals = pd.Series([0,1,1,1,1,1,1,3,3, 3,3,3,3], index=columns)
    df1_new = df[columns].round(decimals)
    df1_new.columns = new_cols
    
    for idx, row in df1_new.iterrows():
        m = row['Method']
        df1_new.at[idx,'Method']=dict_methods[m]
        
    return df1_new

###################################################################################################

# main function to run both evaluations
    
###################################################################################################   

def get_pvalues_single(filename, n=2000, reference_genes=cosmic_genes,pvalue_threshold=0.05,pvalue_position=3, ref_edges=edges,randiter=100):
    """main function to evaluate methods based on their pvalues.
    n: chunk size (genes to consider for each iteration),
    reference_genes: known driver genes (default is CGC),
    pvalue_threhsold: significance threhsold,
    pvalue_position: column number from mutex result file where p-values are stored. Default is 3,
    ref_edges: PPI network edges,
    randiter: amount of iterations to account for randomization.
    """
    genes = get_genes(filename=filename) 
    print('Total Genes:',len(genes))
    
    ## dictionaries for pos/neg
    #1
    
    dict_pairs_cgcg_for_nnb = {}
    dict_nnb_sigLHS_nonsigRHS = {} 
    dict_nnb_sigLHS_sigRHS_LHS = {} 
    dict_nnb_sigLHS_sigRHS_RHS = {} 
    dict_nnb_nonsigLHS_sigRHS = {} 
    dict_nnb_nonsigLHS_nonsigRHS_LHS = {} 
    dict_nnb_nonsigLHS_nonsigRHS_RHS = {} 
    dict_nnb_sum_LHS = {}
    dict_nnb_sum_sig_LHS = {}
    dict_nnb_sum_RHS = {}
    dict_nnb_sum_sig_RHS = {}
    
    dict_norm_nnb_sigLHS_nonsigRHS = {} 
    dict_norm_nnb_sigLHS_sigRHS_LHS = {} 
    dict_norm_nnb_sigLHS_sigRHS_RHS = {} 
    dict_norm_nnb_nonsigLHS_sigRHS = {} 
    dict_norm_nnb_nonsigLHS_nonsigRHS_LHS = {} 
    dict_norm_nnb_nonsigLHS_nonsigRHS_RHS = {} 
    
    #2
    
    dict_pairs_cgcg_for_ncgnb = {}
    dict_ncgnb_sigLHS_nonsigRHS = {} 
    dict_ncgnb_sigLHS_sigRHS_LHS = {} 
    dict_ncgnb_sigLHS_sigRHS_RHS = {} 
    dict_ncgnb_nonsigLHS_sigRHS = {} 
    dict_ncgnb_nonsigLHS_nonsigRHS_LHS = {} 
    dict_ncgnb_nonsigLHS_nonsigRHS_RHS = {} 
    dict_ncgnb_sum_sig_LHS = {}
    dict_ncgnb_sum_LHS = {}
    dict_ncgnb_sum_sig_RHS = {} 
    dict_ncgnb_sum_RHS = {}    
    
    dict_norm_ncgnb_sigLHS_nonsigRHS = {} 
    dict_norm_ncgnb_sigLHS_sigRHS_LHS = {} 
    dict_norm_ncgnb_sigLHS_sigRHS_RHS = {} 
    dict_norm_ncgnb_nonsigLHS_sigRHS = {} 
    dict_norm_ncgnb_nonsigLHS_nonsigRHS_LHS = {}
    dict_norm_ncgnb_nonsigLHS_nonsigRHS_RHS = {}
    
    
    ## dictionaries for cg-ncg neighbor degrees
    dict_neighbors_degree_all = {}
    dict_neighbors_cg = {}
    dict_neighbors_ncg = {}
    
    ## filter intact to contain only these genes
    dict_neighbors = get_neighbors(genes,ref_edges=ref_edges)
    
    ## Cohort specific ref (COSMIC) genes
    cohort_ref_genes = set.intersection(set(reference_genes),set(genes),set(dict_neighbors))
    cg_cg_genes = get_cg_cg_genes(cohort_specific_genes=genes, dict_neighbor=dict_neighbors) 
    print('Cosmic Genes:',len(cohort_ref_genes))
    print('CG-CG:',len(cg_cg_genes))

    #group gene into chunks
    group_of_genes = list(chunks(genes, n=n))
    
    #read file and get min nonzero pvalue
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
    
    ## groupwise computations
    count_g = 0
    
    for group in tqdm(group_of_genes,desc='group'):
        
        dict_temp = {g:{} for g in group}
        
        for line in tqdm(lines):
            line = line.strip().split('\t')
        
            g1 = line[1]
            g2 = line[2]
            
            ## added post memo run
            p_temp = float(line[pvalue_position])
            
            if p_temp==0:
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
                

        count_s = []
        for g in tqdm(dict_temp,desc='genes and neighbors'):
            
            if g in reference_genes and g in dict_neighbors:
                
                ### 1 with non neighbors
                temp_pairs,temp_one, temp_two,temp_three,temp_four,temp_five,\
                temp_six, temp_sumLHS, temp_sumsigLHS, temp_sumRHS,temp_sumsigRHS,\
                temp_one_norm, temp_two_norm,temp_three_norm,temp_four_norm,temp_five_norm, temp_six_norm= get_sig_logpval_counts_cgcg_minus_cgnnb_single(dict_temp[g], neighbor_set=dict_neighbors[g],randiter=randiter)
                if not np.isnan(temp_one):
                    dict_pairs_cgcg_for_nnb[g],dict_nnb_sigLHS_nonsigRHS[g], dict_nnb_sigLHS_sigRHS_LHS[g],dict_nnb_sigLHS_sigRHS_RHS[g],\
                    dict_nnb_nonsigLHS_nonsigRHS_LHS[g], dict_nnb_nonsigLHS_nonsigRHS_RHS[g],dict_nnb_nonsigLHS_sigRHS[g],\
                    dict_nnb_sum_LHS[g],dict_nnb_sum_sig_LHS[g],dict_nnb_sum_RHS[g],\
                    dict_nnb_sum_sig_RHS[g]= temp_pairs,temp_one, temp_two,temp_three,temp_four,temp_five,temp_six, temp_sumLHS,temp_sumsigLHS, temp_sumRHS,temp_sumsigRHS
                
                    dict_norm_nnb_sigLHS_nonsigRHS[g], dict_norm_nnb_sigLHS_sigRHS_LHS[g],dict_norm_nnb_sigLHS_sigRHS_RHS[g],\
                    dict_norm_nnb_nonsigLHS_nonsigRHS_LHS[g], dict_norm_nnb_nonsigLHS_nonsigRHS_RHS[g],\
                    dict_norm_nnb_nonsigLHS_sigRHS[g]= temp_one_norm, temp_two_norm,temp_three_norm,temp_four_norm,temp_five_norm,temp_six_norm
                
                ## 2 with non cosmic neighbors
                temp_pairs,temp_one, temp_two,temp_three,temp_four,temp_five,\
                temp_six, temp_sumLHS, temp_sumsigLHS, temp_sumRHS,temp_sumsigRHS,\
                temp_one_norm, temp_two_norm,temp_three_norm,temp_four_norm,temp_five_norm,\
                temp_six_norm = get_sig_logpval_counts_cgcg_minus_cgncgnb_single(dict_temp[g], neighbor_set=dict_neighbors[g],randiter=randiter)
                if not np.isnan(temp_one):
                
                    dict_pairs_cgcg_for_ncgnb[g],dict_ncgnb_sigLHS_nonsigRHS[g], dict_ncgnb_sigLHS_sigRHS_LHS[g],dict_ncgnb_sigLHS_sigRHS_RHS[g],\
                    dict_ncgnb_nonsigLHS_nonsigRHS_LHS[g], dict_ncgnb_nonsigLHS_nonsigRHS_RHS[g],dict_ncgnb_nonsigLHS_sigRHS[g],\
                    dict_ncgnb_sum_LHS[g],dict_ncgnb_sum_sig_LHS[g],dict_ncgnb_sum_RHS[g],\
                    dict_ncgnb_sum_sig_RHS[g]= temp_pairs,temp_one, temp_two,temp_three,temp_four,temp_five,temp_six, temp_sumLHS, temp_sumsigLHS, temp_sumRHS,temp_sumsigRHS
        
                    dict_norm_ncgnb_sigLHS_nonsigRHS[g], dict_norm_ncgnb_sigLHS_sigRHS_LHS[g],dict_norm_ncgnb_sigLHS_sigRHS_RHS[g],\
                    dict_norm_ncgnb_nonsigLHS_nonsigRHS_LHS[g], dict_norm_ncgnb_nonsigLHS_nonsigRHS_RHS[g],\
                    dict_norm_ncgnb_nonsigLHS_sigRHS[g]= temp_one_norm, temp_two_norm,temp_three_norm,temp_four_norm,temp_five_norm,temp_six_norm
        
                    
        for g in tqdm(dict_neighbors,desc='degrees'):
            dict_neighbors_degree_all[g] = len([v for v in dict_neighbors[g]])
            dict_neighbors_cg[g] = len([v for v in dict_neighbors[g] if v in reference_genes])
            dict_neighbors_ncg[g] = len([v for v in dict_neighbors[g] if v not in reference_genes])
            
    sample_count_2_4 = len(cg_cg_genes)-count_g
    print('sample count',sample_count_2_4)

    return len(cohort_ref_genes), sample_count_2_4,\
    dict_neighbors_degree_all, dict_neighbors_cg, dict_neighbors_ncg, \
    dict_pairs_cgcg_for_nnb,dict_pairs_cgcg_for_ncgnb,\
    dict_nnb_sigLHS_nonsigRHS, dict_nnb_sigLHS_sigRHS_LHS,dict_nnb_sigLHS_sigRHS_RHS,\
    dict_nnb_nonsigLHS_nonsigRHS_LHS, dict_nnb_nonsigLHS_nonsigRHS_RHS,dict_nnb_nonsigLHS_sigRHS,\
    dict_nnb_sum_LHS,dict_nnb_sum_sig_LHS,dict_nnb_sum_RHS, dict_nnb_sum_sig_RHS,\
    dict_norm_nnb_sigLHS_nonsigRHS, dict_norm_nnb_sigLHS_sigRHS_LHS,dict_norm_nnb_sigLHS_sigRHS_RHS,\
    dict_norm_nnb_nonsigLHS_nonsigRHS_LHS, dict_norm_nnb_nonsigLHS_nonsigRHS_RHS,dict_norm_nnb_nonsigLHS_sigRHS,\
    dict_ncgnb_sigLHS_nonsigRHS, dict_ncgnb_sigLHS_sigRHS_LHS,dict_ncgnb_sigLHS_sigRHS_RHS,\
    dict_ncgnb_nonsigLHS_nonsigRHS_LHS, dict_ncgnb_nonsigLHS_nonsigRHS_RHS,dict_ncgnb_nonsigLHS_sigRHS,\
    dict_ncgnb_sum_LHS,dict_ncgnb_sum_sig_LHS,dict_ncgnb_sum_RHS, dict_ncgnb_sum_sig_RHS,\
    dict_norm_ncgnb_sigLHS_nonsigRHS, dict_norm_ncgnb_sigLHS_sigRHS_LHS,dict_norm_ncgnb_sigLHS_sigRHS_RHS,\
    dict_norm_ncgnb_nonsigLHS_nonsigRHS_LHS, dict_norm_ncgnb_nonsigLHS_nonsigRHS_RHS,dict_norm_ncgnb_nonsigLHS_sigRHS    

#Single Neighbor Analysis
# The algorithms are run and the results are stored in dictionaries

list_vals_cgcg_cgnnb = []
list_vals_cgcg_cgnnb_subset = []
list_vals_cgcg_cgncgnb = []

cg_size = {}
cg_cg_size = {} #for getting cg-cg minus cg-ncg size
common_genes = {}
dict_neighbor_degree_all = {}
dict_neighbor_cg={}
dict_neighbor_ncg={}

## 1

dict_pairs_cg_for_cgnnb = {}
dict_cgnnb_sigLHS_nonsigRHS = {} 
dict_cgnnb_sigLHS_sigRHS_LHS = {} 
dict_cgnnb_sigLHS_sigRHS_RHS = {} 
dict_cgnnb_nonsigLHS_sigRHS = {} 
dict_cgnnb_nonsigLHS_nonsigRHS_LHS = {} 
dict_cgnnb_nonsigLHS_nonsigRHS_RHS = {} 
dict_cgnnb_sum_LHS = {}
dict_cgnnb_sumsig_LHS = {}
dict_cgnnb_sum_RHS = {}
dict_cgnnb_sumsig_RHS = {}

dict_norm_cgnnb_sigLHS_nonsigRHS = {} 
dict_norm_cgnnb_sigLHS_sigRHS_LHS = {} 
dict_norm_cgnnb_sigLHS_sigRHS_RHS = {} 
dict_norm_cgnnb_nonsigLHS_sigRHS = {} 
dict_norm_cgnnb_nonsigLHS_nonsigRHS_LHS = {} 
dict_norm_cgnnb_nonsigLHS_nonsigRHS_RHS = {} 

## 2

dict_pairs_cg_for_cgncgnb = {}
dict_cgncgnb_sigLHS_nonsigRHS = {} 
dict_cgncgnb_sigLHS_sigRHS_LHS = {} 
dict_cgncgnb_sigLHS_sigRHS_RHS = {} 
dict_cgncgnb_nonsigLHS_sigRHS = {} 
dict_cgncgnb_nonsigLHS_nonsigRHS_LHS = {} 
dict_cgncgnb_nonsigLHS_nonsigRHS_RHS = {} 
dict_cgncgnb_sum_LHS = {}
dict_cgncgnb_sumsig_LHS = {}
dict_cgncgnb_sum_RHS = {}
dict_cgncgnb_sumsig_RHS = {}

dict_norm_cgncgnb_sigLHS_nonsigRHS = {} 
dict_norm_cgncgnb_sigLHS_sigRHS_LHS = {} 
dict_norm_cgncgnb_sigLHS_sigRHS_RHS = {} 
dict_norm_cgncgnb_nonsigLHS_sigRHS = {} 
dict_norm_cgncgnb_nonsigLHS_nonsigRHS_LHS = {} 
dict_norm_cgncgnb_nonsigLHS_nonsigRHS_RHS = {} 

randiter = 100
randiter_outer = 100
zero_threshold = 0

cols_sum = ['method', 'pairs','case1', 'case2', 'case3', 'case4', 'case5', 'case6', '(1+2)','(3+6)', 'sgm_allCGNB','sgm_sigCGNB','avg_allCGNB','avg_CGNNB','sgm_allCGNNB','sgm_sigCGNNB',\
            'TP', 'FP', 'TN', 'FN', 'PR','SN(TPR)', 'SP(TNR)', 'F1']
cols_sum_ncgnb = ['method', 'pairs','case1', 'case2', 'case3', 'case4', 'case5', 'case6', '(1+2)','(3+6)', 'sgm_allCGNB','sgm_sigCGNB','avg_allCGNB','avg_sigCGNCGB','sgm_allNCGNB','sgm_sigNCGNB',\
                  'TP', 'FP', 'TN', 'FN', 'PR','SN(TPR)', 'SP(TNR)', 'F1']
for m in tqdm(methods):
    print(m)
    filename = dict_infile[m]
    if m=='wext':
        pval_position = 3
    else:
        pval_position=3

    cg_size[m],cg_cg_size[m], dict_neighbor_degree_all[m],dict_neighbor_cg[m], dict_neighbor_ncg[m], \
    dict_pairs_cg_for_cgnnb[m],dict_pairs_cg_for_cgncgnb[m],\
    dict_cgnnb_sigLHS_nonsigRHS[m], dict_cgnnb_sigLHS_sigRHS_LHS[m],dict_cgnnb_sigLHS_sigRHS_RHS[m],\
    dict_cgnnb_nonsigLHS_nonsigRHS_LHS[m], dict_cgnnb_nonsigLHS_nonsigRHS_RHS[m],dict_cgnnb_nonsigLHS_sigRHS[m],\
    dict_cgnnb_sum_LHS[m],dict_cgnnb_sumsig_LHS[m],dict_cgnnb_sum_RHS[m],dict_cgnnb_sumsig_RHS[m],\
    dict_norm_cgnnb_sigLHS_nonsigRHS[m], dict_norm_cgnnb_sigLHS_sigRHS_LHS[m],dict_norm_cgnnb_sigLHS_sigRHS_RHS[m],\
    dict_norm_cgnnb_nonsigLHS_nonsigRHS_LHS[m], dict_norm_cgnnb_nonsigLHS_nonsigRHS_RHS[m],dict_norm_cgnnb_nonsigLHS_sigRHS[m],\
    dict_cgncgnb_sigLHS_nonsigRHS[m], dict_cgncgnb_sigLHS_sigRHS_LHS[m],dict_cgncgnb_sigLHS_sigRHS_RHS[m],\
    dict_cgncgnb_nonsigLHS_nonsigRHS_LHS[m], dict_cgncgnb_nonsigLHS_nonsigRHS_RHS[m],dict_cgncgnb_nonsigLHS_sigRHS[m],\
    dict_cgncgnb_sum_LHS[m],dict_cgncgnb_sumsig_LHS[m],dict_cgncgnb_sum_RHS[m],dict_cgncgnb_sumsig_RHS[m],\
    dict_norm_cgncgnb_sigLHS_nonsigRHS[m], dict_norm_cgncgnb_sigLHS_sigRHS_LHS[m],dict_norm_cgncgnb_sigLHS_sigRHS_RHS[m],\
    dict_norm_cgncgnb_nonsigLHS_nonsigRHS_LHS[m], dict_norm_cgncgnb_nonsigLHS_nonsigRHS_RHS[m],\
    dict_norm_cgncgnb_nonsigLHS_sigRHS[m] = get_pvalues_single(filename,randiter=randiter,n=3000,pvalue_position=pval_position)
    
    print()
    ## for CGNB and CG NNB
    case1_cgnnb = sum(dict_cgnnb_sigLHS_nonsigRHS[m].values())
    case2_cgnnb = sum(dict_cgnnb_sigLHS_sigRHS_LHS[m].values())
    case3_cgnnb = sum(dict_cgnnb_sigLHS_sigRHS_RHS[m].values())
    case4_cgnnb = sum(dict_cgnnb_nonsigLHS_nonsigRHS_LHS[m].values())
    case5_cgnnb = sum(dict_cgnnb_nonsigLHS_nonsigRHS_RHS[m].values())
    case6_cgnnb = sum(dict_cgnnb_nonsigLHS_sigRHS[m].values())
    sumsig_LHS_cgnnb = sum(dict_cgnnb_sumsig_LHS[m].values())
    sum_LHS_cgnnb = sum(dict_cgnnb_sum_LHS[m].values())
    medsig_LHS_cgnnb = np.median(list(dict_cgnnb_sumsig_LHS[m].values())) 
    sumsig_RHS_cgnnb = sum(dict_cgnnb_sumsig_RHS[m].values())
    sum_RHS_cgnnb = sum(dict_cgnnb_sum_RHS[m].values()) 
    medsig_RHS_cgnnb = np.median(list(dict_cgnnb_sumsig_RHS[m].values()))
    
    TP_cgnnb = case1_cgnnb + case2_cgnnb +case3_cgnnb 
    FP_cgnnb = case2_cgnnb + case3_cgnnb +case6_cgnnb
    TN_cgnnb = case1_cgnnb + case4_cgnnb +case5_cgnnb 
    FN_cgnnb = case4_cgnnb + case5_cgnnb +case6_cgnnb
    sensitivity_cgnnb = TP_cgnnb/(TP_cgnnb+FN_cgnnb)
    specificity_cgnnb = TN_cgnnb/(TN_cgnnb+FP_cgnnb)
    precision_cgnnb = TP_cgnnb/(TP_cgnnb+FP_cgnnb)
    f1_score_cgnnb = 2*precision_cgnnb*sensitivity_cgnnb/(precision_cgnnb+sensitivity_cgnnb)

    total_cg_pairs_for_nnb = np.sum(list(dict_pairs_cg_for_cgnnb[m].values()))
    list_vals_cgcg_cgnnb.append([m,total_cg_pairs_for_nnb]
        + [case1_cgnnb, case2_cgnnb,case3_cgnnb, case4_cgnnb, case5_cgnnb, case6_cgnnb,\
           case1_cgnnb + case2_cgnnb, case3_cgnnb + case6_cgnnb,\
           sum_LHS_cgnnb,sumsig_LHS_cgnnb,sum_LHS_cgnnb/total_cg_pairs_for_nnb,sum_RHS_cgnnb/total_cg_pairs_for_nnb,\
           sum_RHS_cgnnb,sumsig_RHS_cgnnb,\
          TP_cgnnb,FP_cgnnb,TN_cgnnb,FN_cgnnb, precision_cgnnb,sensitivity_cgnnb, specificity_cgnnb, f1_score_cgnnb])

    ## for CGNB and NCG NB

    case1_cgncgnb = sum(dict_cgncgnb_sigLHS_nonsigRHS[m].values())
    case2_cgncgnb = sum(dict_cgncgnb_sigLHS_sigRHS_LHS[m].values())
    case3_cgncgnb = sum(dict_cgncgnb_sigLHS_sigRHS_RHS[m].values())
    case4_cgncgnb = sum(dict_cgncgnb_nonsigLHS_nonsigRHS_LHS[m].values())
    case5_cgncgnb = sum(dict_cgncgnb_nonsigLHS_nonsigRHS_RHS[m].values())
    case6_cgncgnb = sum(dict_cgncgnb_nonsigLHS_sigRHS[m].values())
    sumsig_LHS_cgncgnb = sum(dict_cgncgnb_sumsig_LHS[m].values())
    sum_LHS_cgncgnb = sum(dict_cgncgnb_sum_LHS[m].values())
    medsig_LHS_cgncgnb = np.median(list(dict_cgncgnb_sumsig_LHS[m].values()))
    sumsig_RHS_cgncgnb = sum(dict_cgncgnb_sumsig_RHS[m].values())
    sum_RHS_cgncgnb = sum(dict_cgncgnb_sum_RHS[m].values())
    medsig_RHS_cgncgnb = np.median(list(dict_cgncgnb_sumsig_RHS[m].values()))
    
    TP_cgncgnb = case1_cgncgnb + case2_cgncgnb +case3_cgncgnb 
    FP_cgncgnb = case2_cgncgnb + case3_cgncgnb +case6_cgncgnb
    TN_cgncgnb = case1_cgncgnb + case4_cgncgnb +case5_cgncgnb 
    FN_cgncgnb = case4_cgncgnb + case5_cgncgnb +case6_cgncgnb
    sensitivity_cgncgnb = TP_cgncgnb/(TP_cgncgnb+FN_cgncgnb)
    specificity_cgncgnb = TN_cgncgnb/(TN_cgncgnb+FP_cgncgnb)
    precision_cgncgnb = TP_cgncgnb/(TP_cgncgnb+FP_cgncgnb)
    f1_score_cgncgnb = 2*precision_cgncgnb*sensitivity_cgncgnb/(precision_cgncgnb+sensitivity_cgncgnb)
 
    total_cg_pairs_for_ncgnb = np.sum(list(dict_pairs_cg_for_cgncgnb[m].values()))

    list_vals_cgcg_cgncgnb.append([m,total_cg_pairs_for_ncgnb]
        + [case1_cgncgnb, case2_cgncgnb,case3_cgncgnb, case4_cgncgnb, case5_cgncgnb, case6_cgncgnb,\
           case1_cgncgnb + case2_cgncgnb, case3_cgncgnb + case6_cgncgnb,\
           sum_LHS_cgncgnb,sumsig_LHS_cgncgnb,sum_LHS_cgncgnb/total_cg_pairs_for_ncgnb,sum_RHS_cgncgnb/total_cg_pairs_for_ncgnb,\
           sum_RHS_cgncgnb, sumsig_RHS_cgncgnb,\
          TP_cgncgnb,FP_cgncgnb,TN_cgncgnb,FN_cgncgnb,precision_cgncgnb,sensitivity_cgncgnb, specificity_cgncgnb,f1_score_cgncgnb])
    
df_summary_cgcg_cgnnb = pd.DataFrame(list_vals_cgcg_cgnnb, columns = cols_sum)
df_summary_cgcg_cgncgnb = pd.DataFrame(list_vals_cgcg_cgncgnb, columns = cols_sum_ncgnb)

df_summary_cgcg_cgnnb=transaction(df_summary_cgcg_cgnnb)
df_summary_cgcg_cgncgnb=transaction(df_summary_cgcg_cgncgnb)


if not os.path.exists(save_path+ '/results_counts_eval1_cgcg_cgnnb_tpfp/'):
    os.makedirs(save_path+'/results_counts_eval1_cgcg_cgnnb_tpfp/')
if not os.path.exists(save_path+ '/results_counts_eval2_cgcg_cgncgnb_tpfp/'):
    os.makedirs(save_path+'/results_counts_eval2_cgcg_cgncgnb_tpfp/')

outfile_cgcg_cgnnb = save_path+'/results_counts_eval1_cgcg_cgnnb_tpfp/{}_t{}_{}.txt'.format(c,t, '_'.join(methods))
outfile_cgcg_cgncgnb = save_path+'/results_counts_eval2_cgcg_cgncgnb_tpfp/{}_t{}_{}.txt'.format(c,t, '_'.join(methods))
df_summary_cgcg_cgnnb.to_csv(outfile_cgcg_cgnnb, index=False, sep='\t')
df_summary_cgcg_cgncgnb.to_csv(outfile_cgcg_cgncgnb, index=False, sep='\t')

