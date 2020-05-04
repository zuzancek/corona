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
import init_common as x
import helpers as hp
from random import sample
from sklearn.utils import shuffle

sns.set(rc={'figure.figsize':(11, 4)})

R0_target = 2.2
N_nodes = 1000000

## 1. add distribution for recovery time / gamma
Trec_mean  = 3.8 #6.5 is serial interval, this is very conservative estimate of infectious period, starts approx 2.5 days before the end of inc.period
Trec_sig = 0.62 
Tinf_mean = 3.8
m0 = Trec_mean*(Trec_sig**2)
s0 = 1/(Trec_sig**2)
Trec_vec = np.random.gamma(m0,s0,N_nodes)
Tinf_vec = np.random.gamma(Tinf_mean/s0,s0,N_nodes)
Trec_vec_list = Trec_vec.tolist()
Tinf_vec_list = Tinf_vec.tolist()
Trec_mean = np.average(Trec_vec)
gamma_mean = 1/Trec_mean 
gamma_vec = 1/Trec_vec
## use this -->>
gamma_vec_list = gamma_vec.tolist()


## 2.add distribution for R0 
## 2.A self-isolated
isol_share = 0.7
N_nodes_isol = round(isol_share*N_nodes)
R0_scale = 100 #10**2
R0_isol_mean = 0.2*3.96# 1.43
R0_isol_vec = np.random.gamma(R0_scale*R0_isol_mean,1/R0_scale,N_nodes_isol)

## 2.B non-isolated
a0 = R0_isol_mean+0*np.std(R0_isol_vec)#2
a1 = 10
scale_nonisol = 2.75#1.65
N_nodes_nonisol = N_nodes-N_nodes_isol
R0_nonisol_vec = hp.power_law_pdf(np.random.uniform(0,1,N_nodes_nonisol),a0,a1,scale_nonisol)

## 2.C joint distribution
R0_vec = shuffle(np.transpose([*np.transpose(R0_isol_vec),*np.transpose(R0_nonisol_vec)]),random_state=0)
R0_vec_list = R0_vec.tolist()
R0_mean = np.mean(R0_vec)
R0_std = np.std(R0_vec)

## 3. exctract Beta
beta_vec = R0_vec/Trec_vec
beta_vec_mean = np.mean(beta_vec)
## use this -->>
beta_vec_list = beta_vec.tolist()

with open('./src/OD.pickle','rb') as f:
    OD=pickle.load(f)
np.fill_diagonal(OD,0)
for j in range(x.N_locs):
    kappa = np.sum(OD[:,j])/x.N_popul[j]
    if kappa>1:
        OD[:,j] = OD[:,j]/kappa


def get_vectors(R0_scale=1):
    beta_vec_new = R0_scale*beta_vec
    beta_vec_list_new = beta_vec_new.tolist()
    beta_0 = R0_scale*np.random.gamma(R0_scale*R0_isol_mean,1/R0_scale,N_nodes)/Trec_vec
    beta_0_list = beta_0.tolist()
    beta_1 = R0_scale*hp.power_law_pdf(np.random.uniform(0,1,N_nodes),a0,a1,scale_nonisol)/Trec_vec
    beta_1_list = beta_vec_list_new
    return beta_vec_list_new,Trec_vec_list,beta_0_list,beta_1_list

