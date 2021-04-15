# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 20:26:00 2021

@author: cansu yalçın
"""
import zipfile
with zipfile.ZipFile('data/fishers_mutation_filtered_ep_data.zip', 'r') as zip_ref:
    zip_ref.extractall("data/")
    
with zipfile.ZipFile('data/megsa_mutation_filtered_ep_data.zip', 'r') as zip_ref:
    zip_ref.extractall("data/")

with zipfile.ZipFile('data/memo_mutation_filtered_ep_data.zip', 'r') as zip_ref:
    zip_ref.extractall("data/")
    
with zipfile.ZipFile('data/wext_mutation_filtered_ep_data.zip', 'r') as zip_ref:
    zip_ref.extractall("data/")
    
with zipfile.ZipFile('data/discover_mutation_filtered_ep_data/discover_mutation_filtered_ep_data_normal.zip', 'r') as zip_ref:
    zip_ref.extractall("data/discover_mutation_filtered_ep_data")

with zipfile.ZipFile('data/discover_mutation_filtered_ep_data/discover_mutation_filtered_ep_data_stratified.zip', 'r') as zip_ref:
    zip_ref.extractall("data/discover_mutation_filtered_ep_data")