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
import init as x

def power_law_pdf(y,a0,a1,k):
    z = ((a1**(1.0-k)-a0**(1.0-k))*y+a0**(1.0-k))**(1.0/(1.0-k))
    return z

def adjust_beta(beta1,beta0,do_print=False):
    beta0list = np.array(beta0)
    beta1 = (beta1/np.average(beta0list))*beta0list
    if do_print:
        print('new beta mean = ',np.average(beta1))
    return beta1.to_list()

## Funkcia pre výpočet priemeru a odchylky zo simulácií
def mean_list(x,colname):
    z = pd.DataFrame(x[0][colname])
    if len(z)>1:
        for idx in x[1:]:
            z = pd.concat([z,pd.DataFrame(idx[colname])],1)
    return z.apply(np.mean,1)

def std_list(x,colname):
    z = pd.DataFrame(x[0][colname])
    if len(z)>1:
        for idx in x[1:]:
            z = pd.concat([z,pd.DataFrame(idx[colname])],1)
    return z.apply(np.std,1)

def mean_diff_list(x,tidx):
    z = pd.DataFrame(x[0][:,tidx,:].sum(0)[1:]-x[0][:,tidx,:].sum(0)[:-1])
    if len(x)>1:
        for idx in x[1:]:
              z = pd.concat([z,pd.DataFrame(idx[:,tidx,:].sum(0)[1:]-idx[:,tidx,:].sum(0)[:,-1])],1)
    return z.apply(np.mean,1)

def std_diff_list(x,tidx):
    z = pd.DataFrame(x[0][:,tidx,:].sum(0)[1:]-x[0][:,tidx,:].sum(0)[:-1])
    if len(x)>1:
        for idx in x[1:]:
              z = pd.concat([z,pd.DataFrame(idx[:,tidx,:].sum(0)[1:]-idx[:,tidx,:].sum(0)[:,-1])],1)
    return z.apply(np.std,1)

# OD matrix
# load data and process matrix
def get_OD_matrix():
    with open('./src/OD.pickle','rb') as f:
        OD=pickle.load(f)
        f.close()
    #np.fill_diagonal(OD,0)
    #for idx in range(N_locs):
    #    sh = np.sum(OD[:,idx])/N_popul[idx]
    #    if sh>1:
    #        OD[:,idx] = OD[:,idx]/sh
    return OD


def setup_paths(R0_type,alpha_):
    
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
    if R0_type == 1:
        out_filename_ext = out_filename_ext+"R0low"
        out_fig_ext = out_fig_ext+"R0low"
        out_stat_ext = out_stat_ext+"R0low"
    elif R0_type == 2:
        out_filename_ext = out_filename_ext+"R0high"
        out_fig_ext = out_fig_ext+"R0high"
        out_stat_ext = out_stat_ext+"R0high"
        
    out_filename_root = out_filename_root+"/"+out_filename_ext
    out_fig_root = out_fig_root+"/"+out_fig_ext
    out_stat_root = out_stat_root+"/"+out_stat_ext
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
    return out_filename_root,out_fig_root,out_stat_root  
    
    

