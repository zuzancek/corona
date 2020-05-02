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
np.random.seed(1000)

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
# corrrection for unobserved    
first_infections_original=first_infections
first_infections_correction_multiplier = 16
first_infections=first_infections_original*first_infections_correction_multiplier/4
unobs_infection_num = int((first_infections_correction_multiplier*3/4)*sum(first_infections_original))
unobs_idx = np.random.randint(0,N_locs,unobs_infection_num)
for i in np.arange(unobs_infection_num):
    ri = int(np.random.uniform(1)*N_locs)
    first_infections[unobs_idx[i]]=first_infections[unobs_idx[i]]+1
tau = 16/24
R0_default = 1.48#1.913 #1.47

## technical params
N_per = 200  # number of simulated periods (days)  
N_simul = 64  # number of repetitions/independent runs

global fnc_type
global R0_type
fnc_type = 0
R0_type = 0
# OD = get_OD_matrix()

with open('./src/OD_old.pickle','rb') as f:
    OD=pickle.load(f)

## paths and directories
out_filename_root = "./out"
out_fig_root = "./fig"
out_stat_root = "./stat"
out_filename = "simul_SIR"
out_filename_raw0 = "simul_SIR_raw"

## epidemiology parameters
N_nodes = 1000000
Tinc_mean = 5.1
T_sig = 0.62
Tinc_vec = np.random.gamma(Tinc_mean*(T_sig**2),1/(T_sig**2),N_nodes)
Tinc_inv_vec = 1/Tinc_vec
Tinc_inv_list = Tinc_inv_vec.tolist()
SI_mean = 6.5
presymp_period = 0.5
Tinf_mean = SI-(Tinc_mean-presymp_period)
Tinf_vec = np.random.gamma(Tinf_mean*(T_sig**2),1/(T_sig**2),N_nodes)
Tinf_inv_vec = 1/Tinf_vec
Tinf_inv_list = Tinf_inv_vec.tolist()
## distribution for R0 
R0_target = 2.2
# self-isolated
isol_share = 0.7
N_nodes_isol = round(isol_share*N_nodes)
R0_scale = 100 #10**2
R0_isol_mean = 0.25*3.96
R0_isol_vec = np.random.gamma(R0_scale*R0_isol_mean,1/R0_scale,N_nodes_isol)
# non-isolated
a0 = R0_isol_mean+np.std(R0_isol_vec)#2
a1 = 14#16
scale_nonisol = 2.3#1.65
N_nodes_nonisol = N_nodes-N_nodes_isol
R0_nonisol_vec = hp.power_law_pdf(np.random.uniform(0,1,N_nodes_nonisol),a0,a1,scale_nonisol)
# joint distribution
R0_vec = shuffle(np.transpose([*np.transpose(R0_isol_vec),*np.transpose(R0_nonisol_vec)]),random_state=0)
R0_vec_list = R0_vec.tolist()
R0_mean = np.mean(R0_vec)
R0_std = np.std(R0_vec)
# exctract Beta
beta_vec = R0_vec/Trec_vec
beta_vec_mean = np.mean(beta_vec)
beta_list = beta_vec.tolist()

## policy params
rho_trans = 1 # default

## clinical parameters:
# probabilities & shares
omega_obs = 1.0/first_infections_correction_multiplier
omega_asymp = .5
prob_hosp = .2
prob_icu = .05/prob_hosp
prob_die = 0.02/(prob_icu*prob_hosp)
# times
T_hosp = 5
Trec_hosp = 10
T_icu=3
Trec_icu = 10
T_death = 8
Trec_mild = 8    