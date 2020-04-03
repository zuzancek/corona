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
Trec_mean  = 6.5 # serial interval
Trec_sig = 0.62
m0 = Trec_mean*(Trec_sig**2)
s0 = 1/(Trec_sig**2)
Trec_vec = np.random.gamma(m0,s0,N_nodes)
Trec_mean = np.average(Trec_vec)
gamma_mean = 1/Trec_mean 


## 2.add distribution for beta 
## 2.A self-isolated
isol_share = 0.7
N_nodes_isol = round(isol_share*N_nodes)
# their beta follows gamma dist. with very small variance
beta_scale = 1000
beta_isol_mean = 0.15
beta_isol = beta_scale*beta_isol_mean
beta_vec_isol = np.random.gamma(beta_isol,1.0/beta_scale,N_nodes_isol)
Trec_isol = Trec_vec[0:N_nodes_isol]
R0_isol_vec = beta_vec_isol*Trec_isol

## 2.B non-isolated
a0 = 0.15
a1 = 10
scale_nonisol = 2.45
N_nodes_nonisol = N_nodes-N_nodes_isol
beta_vec_nonisol = hp.power_law_pdf(np.random.uniform(0,1,N_nodes_nonisol),a0,a1,scale_nonisol)

Trec_nonisol = Trec_vec[N_nodes_isol:N_nodes]
R0_nonisol_vec = beta_vec_nonisol*Trec_nonisol

## 3. joint distribution
beta_vec = shuffle(np.transpose([*np.transpose(beta_vec_isol),*np.transpose(beta_vec_nonisol)]),random_state=0)
beta_vec_list = beta_vec.tolist()

R0_vec = beta_vec*Trec_vec
R0_mean = np.mean(R0_vec)
R0_std = np.std(R0_vec)
