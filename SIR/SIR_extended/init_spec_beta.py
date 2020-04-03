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
# Trec is deterministic in fact
Trec_mean  = 6.5 
Trec_vec = np.full(N_nodes,Trec_mean)
gamma_mean = 1/Trec_mean 
gamma_vec = 1/rec_vec
## use this -->>
gamma_vec_list = gamma_vec.tolist()


## 2.add distribution for R0 
## 2.A self-isolated
isol_share = 0.7
N_nodes_isol = round(isol_share*N_nodes)
R0_scale = 64 #10**2
R0_isol_mean = .97
R0_isol_vec = np.random.gamma(R0_scale*R0_isol_mean,1/R0_scale,N_nodes_isol)

## 2.B non-isolated
a0 = R0_isol_mean+2*np.std(R0_isol_vec)
a1 = 14
scale_nonisol = 2
N_nodes_nonisol = N_nodes-N_nodes_isol
R0_nonisol_vec = hp.power_law_pdf(np.random.uniform(0,1,N_nodes_nonisol),a0,a1,scale_nonisol)

## 2.C joint distribution
R0_vec = shuffle(np.transpose([*np.transpose(R0_isol_vec),*np.transpose(R0_nonisol_vec)]),random_state=0)
R0_vec_list = R0_vec.tolist()
R0_mean = np.mean(R0_vec)
R0_std = np.std(R0_vec)

## 3. exctract Beta
Beta_vec = R0_vec/Trec_vec
Beta_vec = Beta_vec[Beta_vec<2]
## use this -->>
beta_vec_list = beta.tolist()

