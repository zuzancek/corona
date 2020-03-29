import sys
import numpy as np
from tqdm import tqdm_notebook
import pickle
import pandas as pd
import plotly
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
import time
import os
sns.set(rc={'figure.figsize':(11, 4)})

## Funkcia pre výpočet priemeru zo simulácií
def sumlist(x):
    tmp=x[0]
    for i in x[1:]:
        tmp=tmp+i
    return tmp/len(x)

## Údaje k počtu obyvateľov na obec 
pop = pd.read_excel('./src/munic_pop.xlsx')
pop_N = np.array(pop['popul'])
N_popul = pop.popul.to_numpy()        # Populacia (vektor)
N_popul_size = np.sum(N_popul)         # pocet obyvatelov
N_locs = len(N_popul)                 # Pocet obci

## Priradenie GPS suradnic pre kazdu obec
def get_coors_long(x):
    return float(df_coords.loc[df_coords.IDN4.apply(str)==x,'long'])

def get_coors_lat(x):
    return float(df_coords.loc[df_coords.IDN4.apply(str)==x,'lat'])
df_coords=pd.read_excel('./src/obce1.xlsx')
data_i=pop
data_i.loc[:,'long']=data_i.munic.apply(str).apply(get_coors_long)
data_i.loc[:,'lat']=data_i.munic.apply(str).apply(get_coors_lat)

# OD matrix
# load data and process matrix
def get_OD_matrix():
    with open('./src/OD.pickle','rb') as f:
        OD=pickle.load(f)
        f.close()
    np.fill_diagonal(OD,0)
    for idx in range(N_locs):
        sh = np.sum(OD[:,idx])/N_popul[idx]
        if sh>1:
            OD[:,idx] /= sh
    return OD

## load first infections
nakazy_sk = pd.read_excel('./src/cases.xlsx')
first_infections=np.zeros(N_locs)
for i in np.arange(nakazy_sk.shape[0]):
    first_infections[pop.munic==nakazy_sk.KOD.iloc[i]]=nakazy_sk.ID.iloc[i]
    
first_infections_original=first_infections
first_infections_correction_multiplier = 6
first_infections=first_infections_original*first_infections_correction_multiplier

## ALTERNATIVE
global R0
R0 = 2.4
Tinf = 3
Tinc = 5
gamma = 1/5#(Tinf+0*Tinc)
global beta
beta = R0*gamma

## technical params
N_k = pop.popul.to_numpy()          # Populacia
locs_len = len(N_k)                 # Pocet obci
simul_len = 200   
simul_cnt = 64
public_trans_high = 1
public_trans_mid = 0.8
public_trans_low = 0.6

global fnc_type
global R0_type
fnc_type = 0
R0_type = 0
OD = get_OD_matrix()

data_senior=pd.read_excel('./src/senior.xlsx')
data_senior.loc[:,'munic']=data_senior.munic.apply(lambda x: x[-6:]).apply(int)
data_senior=data_senior.sort_values(by=['munic'])

N_k_s = N_k-data_senior.senior.to_numpy()
locs_len_s = len(N_k_s)

## paths and directories
out_filename_root = "./out"
out_fig_root = "./fig"
out_stat_root = "./stat"
out_filename = "SIR.pickle"

def setup_paths(fnc_type,R0_type):
    out_filename_root = "./out"
    out_fig_root = "./fig"
    out_stat_root = "./stat"
    try:
        os.mkdir(out_filename_root)
    except:  
        ethrown=True
    try:
        os.mkdir(out_fig_root)
    except:
        ethrown=True
    try:
        os.mkdir(out_stat_root)
    except:
        ethrown=True
    out_filename_ext = ""
    out_fig_ext = ""
    out_stat_ext = ""
    if fnc_type == 0:
        out_filename_dir = ""
        out_fig_dir = ""
        out_stat_dir = ""
        if R0_type == 0:
            return out_filename_root,out_fig_root,out_stat_root
    else:        
        out_filename_ext = "sen"
        out_fig_ext = "sen"
        out_stat_ext = "sen"
    if R0_type == 1:
        out_filename_ext = out_filename_ext+"R0low"
    elif R0_type == 2:
        out_filename_ext = out_filename_ext+"R0high"
    
    out_filename_root = out_filename_root+"/"+out_filename_ext
    try:
        os.mkdir(out_filename_root)
    except:
        ethrown=True
    out_fig_root = out_fig_root+"/"+out_fig_ext
    try:
        os.mkdir(out_fig_root)
    except:
        ethrown=True
    out_stat_root = out_stat_root+"/"+out_stat_ext
    try:
        os.mkdir(out_stat_root)
    except:
        ethrown=True
    return out_filename_root,out_fig_root,out_stat_root  
    
    

