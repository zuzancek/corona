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

# Trec ~ Gamma, Beta ~ Gamma

N_nodes = 10000000
Trec = x.Tinc+x.Tinf
Trec_mean  = 6.5 # serial interval
Trec_sig = 0.62
m0 = Trec_mean*(Trec_sig**2)
s0 = 1/(Trec_sig**2)
rec_vec = np.random.gamma(m0,s0,N_nodes)
Trec_vec_list = rec_vec.tolist()
Trec_mean = np.average(rec_vec)
gamma_mean = 1/Trec_mean 
gamma_vec = 1/rec_vec
gamma_vec_list = gamma_vec.tolist()

beta_mean = x.R0_default*gamma_mean
beta_scale = 1/100
beta_vec = np.random.gamma(beta_mean/beta_scale,beta_scale,N_nodes)
beta_vec_list = beta_vec.tolist()

def get_vectors(R0_scale=1):
    R0 = x.R0_default*R0_scale
    bm = R0*gamma_mean
    bv = np.random.gamma(bm/beta_scale,beta_scale,N_nodes)
    bv_list = bv.tolist()
    return bv_list,Trec_vec_list 

