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

## load first infections
nakazy_sk = pd.read_excel('./src/cases.xlsx')
first_infections=np.zeros(N_locs)
for i in np.arange(nakazy_sk.shape[0]):
    first_infections[pop.munic==nakazy_sk.KOD.iloc[i]]=nakazy_sk.ID.iloc[i]
    
first_infections_original=first_infections
first_infections_correction_multiplier = 6
first_infections=first_infections_original*first_infections_correction_multiplier

# global R0
R0 = 2.2
Tinf = 3
Tinc = 5
gamma = 1/(Tinc+Tinf)

## technical params
N_per = 200  # number of simulated periods (days)  
N_simul = 60 # number of repetitions/independent runs
public_trans_high = 1
public_trans_mid = 0.8
public_trans_low = 0.6

global fnc_type
global R0_type
fnc_type = 0
R0_type = 0
# OD = get_OD_matrix()
with open('./src/OD_old.pickle','rb') as f:
    OD=pickle.load(f)
    f.close()

data_senior=pd.read_excel('./src/senior.xlsx')
data_senior.loc[:,'munic']=data_senior.munic.apply(lambda x: x[-6:]).apply(int)
data_senior=data_senior.sort_values(by=['munic'])

N_popul_nosen = N_popul-data_senior.senior.to_numpy()
N_popul_nosen_size = np.sum(N_popul_nosen)         # pocet obyvatelov bez seniorov
N_locs_nosen = len(N_popul_nosen)

## paths and directories
out_filename_root = "./out"
out_fig_root = "./fig"
out_stat_root = "./stat"
out_filename_pct0 = "simul_SIR_pct"
out_filename_raw0 = "simul_SIR_raw"
