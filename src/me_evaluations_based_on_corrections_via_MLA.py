# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 22:44:42 2021

@author: Rafsan Ahmed & Cansu Yalcin

ME Evaluations Based on Corrections via MLA

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
import os
import gc
from scipy.stats.stats import pearsonr
from itertools import combinations
import random
import multiprocessing
import concurrent.futures
import string
random.seed(1234)

c='COADREAD' # cancer type
t=20 #mutation threshold
methods = ['discover','discover_strat','fishers','megsa','memo','wext']#'discover_strat','discover','discover_strat']#
inpath_mla = '../data/'

#Load MLA
MLA_infile = inpath_mla + 'MLA_ep_mutation_filtered_all_genes/{}_MLA_standardized.txt'.format(c)

#intact 
inpath_intact = '../data/'
intact_edge_file = inpath_intact+'intact_nodupl_edge_file.txt'
intact_index_file = inpath_intact+'intact_nodupl_index_file.txt'

#cosmic
inpath_cosmic = '../data/'
cosmic_infile = inpath_cosmic+'Census_allFri_Apr_26_12_49_57_2019.tsv'

with open(intact_index_file, 'r') as f:
    indices = {line.split()[0]:line.split()[1] for line in f.readlines()}

with open(intact_edge_file, 'r') as f:
    edges = [(indices[line.split()[0]].upper(),indices[line.split()[1]].upper()) for line in f.readlines()]

with open(cosmic_infile,'r') as f:
    cosmic_genes = [line.split()[0].upper() for line in f.readlines()[1:]]
    
with open(MLA_infile, 'r') as f:
    MLA = {line.split()[0]: float(line.split()[1]) for line in f.readlines()}
    
#FUNCTIONS
def set_input_file_paths(methods, c, t):
    '''set input paths as per the location of the MEX results
    m:methods, c: cancer type, t: threshold'''
    
    strat_list=['COADREAD']
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
    
    return dict_infile, dict_infile_intact

def load_cgcg_pairs(methods,dict_infile, ref_genes=cosmic_genes):
    '''
    load mutual exclusivity values for CGC-CGC pairs
    methods: mutual exclusivity methods, ref_genes=reference driver genes (COSMIC)
    '''
    d_out ={}
    for m in tqdm(methods):
        d_nb={}
        infile=dict_infile[m]
        with open(infile) as f: 
            for line in tqdm(f.readlines()): 
                line=line.split()
                g1=line[1] 
                g2=line[2] 
                if g1 in ref_genes and g2 in ref_genes: 
                    if g1 not in d_nb: 
                        d_nb[g1]={} 
                    if g2 not in d_nb:
                        d_nb[g2]={}
                    d_nb[g1][g2]=float(line[3])
                    d_nb[g2][g1]=float(line[3])
        d_out[m]=d_nb

    return d_out

def get_sig_dict(d,m):
    '''get percentage significance from a dictionary of CGC-CGC neighbors
    '''
    d_out={}
    for g in d:
        count=0
        for k,v in d[g].items():
            if m=='wext' and v==0:
                continue
            elif v<0.05:
                count+=1
        d_out[g]= float(count)/float(len(d[g]))
    return d_out
                
def get_sig_dict_from_random_sampling(d,d_nb,method):
    '''get percentage significance through random sampling
    '''
    common_genes = set.intersection(set(d), set(d_nb))
    dict_sig = {}
    for g in common_genes:
        l_temp=[]
        for i in range(100):
            k_temp = {k:d[g][k] for k in random.sample(list(d[g]),k=len(d_nb[g]))}
            
            count=0
            for k,v in k_temp.items():
                if method=='wext' and v==0:
                    continue
                elif v<0.05:
                    count+=1
            
            l_temp.append(float(count)/len(k_temp))
        
        dict_sig[g] = np.mean(l_temp)
    
    return dict_sig

def plot_from_dict(d, d_mla,ax, annotation_threshold_y=-0.0,annotation_threhsold_x=1.5,title=''):
    l_mla_all = []
    l_sig_all = []
    
    l_mla_nonrare = []
    l_sig_nonrare = []
    
    l_mla_rare = []
    l_sig_rare = []
    
    for k,v in d.items():
        if k in d_mla:
            l_mla_all.append(d_mla[k])
            l_sig_all.append(v)
            
                
    ax.scatter(l_mla_all,l_sig_all, c='K',alpha=0.3)
    
    r,p = pearsonr(l_mla_all,l_sig_all)
    print('r:{}\nP:{}'.format(r,p))
    ax.text(0.85, 0.97, 'r = {}\nP = {:.3g}'.format(round(r,2), p), ha='left', va='top', fontsize=8,transform=ax.transAxes)
    
    for k,v in d.items():
        if v>annotation_threshold_y and d_mla[k]<annotation_threhsold_x:
            ax.annotate(k, (d_mla[k],v), ha='left',va='top',xytext=(d_mla[k]+0.08,v-0.01))
    
#     ax.set_xlim(min(l_mla_all)-0.5,max(l_mla_all)+0.5)
    ax.set_xlim(min_mla-0.5,max_mla+0.5)
    ax.set_ylim(-0.05,1.05)        
    labels = [str(int(round(float(item)*100))) for item in ax.get_yticks()]
    print(labels)
    
    ax.set_yticklabels(labels)#
    ax.set_xlabel('Mutation Load Association (MLA)')
    ax.set_ylabel('Percentage of Significant Findings (P<0.05)')
    ax.set_title(title)
    
def plot_from_dict_arrow(d,d_nb, d_mla,ax, annotation_threshold_y=-0.0,annotation_threhsold_x=1.5,title='',c1='C3',c2='K'):
    '''plot with direction of arrow going from dict d to dict d_nb
    annotation thresholds define the points until which gene annotation 
    is done in the graph 
    '''
    l_mla_all = []
    l_sig_all = []
    l_sig_nb = []
    for k,v in d.items():
        if k in d_mla and k in d_nb:
            l_mla_all.append(d_mla[k])
            l_sig_all.append(v)
            l_sig_nb.append(d_nb[k]) 
    ax.scatter(l_mla_all,l_sig_nb, c=c1,alpha=0.7)
    ax.scatter(l_mla_all,l_sig_all, c=c2,alpha=0.2)
    
    r,p = pearsonr(l_mla_all,l_sig_nb)
    print('r:{}\nP:{}'.format(r,p))
    ax.text(0.85, 0.97, 'r = {}\nP = {:.3g}'.format(round(r,2), p), ha='left', va='top', fontsize=8,transform=ax.transAxes)
    
    for k,v in d_nb.items():
        if v>annotation_threshold_y and d_mla[k]<annotation_threhsold_x:
            ax.annotate(k, (d_mla[k],v), ha='left',va='top',xytext=(d_mla[k]+0.08,v-0.01))
            
    for _mla,_all,_nb in zip(l_mla_all, l_sig_all,l_sig_nb):
        if _mla<annotation_threhsold_x and abs(_all-_nb)>0.01:
            if _all<_nb:
                ax.annotate("", xytext=(_mla, _all), xy=(_mla, _all+_nb-_all), arrowprops=dict(arrowstyle="->", alpha=0.2)) 
            else:
                ax.annotate("", xytext=(_mla, _all), xy=(_mla, _all+_nb-_all), arrowprops=dict(arrowstyle="->", alpha=0.2))

    ax.set_xlim(min_mla-0.5,max_mla+0.5)
    ax.set_ylim(-0.05,1.05)        
    labels = [str(int(round(float(item)*100))) for item in ax.get_yticks()]

    ax.set_yticklabels(labels)
    ax.set_xlabel('Mutation Load Association (MLA)')
    ax.set_ylabel('Percentage of Significant Findings (P<0.05)')
    ax.set_title(title)

def get_mla_limits(d_y,d_x):
    ''' get max and min MLA for same x axis
        d_y: sig values
        d_x: MLA values
    '''
    l = []
    for k,v in d_y.items():
        if k in d_x:
            l.append(d_x[k])
    
    return max(l), min(l)

    
dict_infile, dict_infile_intact = set_input_file_paths(methods, c,t)
dict_cg_cg = load_cgcg_pairs(methods,dict_infile)
dict_cg_cg_nb = load_cgcg_pairs(methods,dict_infile_intact)

            
## plot percentage significance

dict_methods = {'discover': 'DISCOVER',
               'discover_strat': 'DISCOVER Strat',
               'fishers': 'Fisher\'s Exact Test',
               'megsa': 'MEGSA',
               'memo': 'MEMO',
               'wext': 'WExT'}


sfx = '../results_main/evaluation_results/percent_sig_figures/'
outpath_cgcg_nb = sfx

if not os.path.exists(outpath_cgcg_nb):
    os.makedirs(outpath_cgcg_nb)

val_x = 2.0
val_y = 0.5
fig, axes = plt.subplots(len(methods),3,figsize=(25,40))
letters = list(string.ascii_uppercase)[:len(methods)*3] 
axiter=axes.flat 
for i,m in tqdm(enumerate(methods)):
    
    if m=='discover_strat':
        dict_sig = get_sig_dict(dict_cg_cg[m],m)
        dict_sig_discover = get_sig_dict(dict_cg_cg['discover'],'discover')
        dict_sig_rand = get_sig_dict_from_random_sampling(dict_cg_cg[m],dict_cg_cg_nb[m],m) 
        dict_sig_rand_discover = get_sig_dict_from_random_sampling(dict_cg_cg['discover'],dict_cg_cg_nb['discover'],'discover')
        dict_sig_cg_cg_nb = get_sig_dict(dict_cg_cg_nb[m], m)

        #main plot
        max_mla, min_mla = get_mla_limits(dict_sig,MLA)
        ax=axiter[3*i]
        letter=letters[3*i]
        title = '{} CGC-CGC pairs'.format(dict_methods[m])
        plot_from_dict_arrow(dict_sig_discover,dict_sig,MLA,ax,title=title,c1='C0', annotation_threhsold_x=val_x, annotation_threshold_y=val_y)

        ax.text(-0.06, 1.1, letter, transform=ax.transAxes,
              fontsize=16, fontweight='bold', va='top', ha='right')
        
        #random sampled plot
        ax=axiter[3*i+1]
        letter=letters[3*i+1]
        title = '{} CGC-CGC pairs (random)'.format(dict_methods[m])
        plot_from_dict_arrow(dict_sig_rand_discover,dict_sig_rand,MLA,ax,title=title,c1='C0', annotation_threhsold_x=val_x, annotation_threshold_y=val_y)

        ax.text(-0.06, 1.1, letter, transform=ax.transAxes,
              fontsize=16, fontweight='bold', va='top', ha='right')
        
        #neighbors plot
        ax=axiter[3*i+2]
        letter=letters[3*i+2]
        title = '{} CGC-CGC neighbors'.format(dict_methods[m])
        plot_from_dict_arrow(dict_sig_rand,dict_sig_cg_cg_nb,MLA,ax,title=title,c2='C0', annotation_threhsold_x=val_x, annotation_threshold_y=val_y)

        ax.text(-0.06, 1.1, letter, transform=ax.transAxes,
              fontsize=16, fontweight='bold', va='top', ha='right')
    
    else:
        dict_sig = get_sig_dict(dict_cg_cg[m], m)
        dict_sig_rand = get_sig_dict_from_random_sampling(dict_cg_cg[m],dict_cg_cg_nb[m],m)
        dict_sig_cg_cg_nb = get_sig_dict(dict_cg_cg_nb[m], m)

        #main plot
        max_mla, min_mla = get_mla_limits(dict_sig,MLA)
        ax=axiter[3*i]
        letter=letters[3*i]
        title = '{} CGC-CGC pairs'.format(dict_methods[m])
        plot_from_dict(dict_sig,MLA,ax,title=title,annotation_threhsold_x=val_x,annotation_threshold_y=val_y)

        ax.text(-0.06, 1.1, letter, transform=ax.transAxes,
              fontsize=16, fontweight='bold', va='top', ha='right')
        
        #randomized plot
        ax=axiter[3*i+1]
        letter=letters[3*i+1]
        title = '{} CGC-CGC pairs (random)'.format(dict_methods[m])
        plot_from_dict(dict_sig_rand,MLA,ax,title=title,annotation_threhsold_x=val_x,annotation_threshold_y=val_y)
        ax.text(-0.06, 1.1, letter, transform=ax.transAxes,
              fontsize=16, fontweight='bold', va='top', ha='right')
         
        ax=axiter[3*i+2]
        letter=letters[3*i+2]
        title = '{} CGC-CGC neighbors'.format(dict_methods[m])
        plot_from_dict_arrow(dict_sig_rand,dict_sig_cg_cg_nb,MLA,ax,title=title, annotation_threhsold_x=val_x, annotation_threshold_y=val_y)

        ax.text(-0.06, 1.1, letter, transform=ax.transAxes,
              fontsize=16, fontweight='bold', va='top', ha='right')

plt.savefig(outpath_cgcg_nb+'{}_t{}_percsig_random_fig.pdf'.format(c,t),format='pdf', bbox_inches='tight')
plt.show()