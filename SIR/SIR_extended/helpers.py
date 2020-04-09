import sys
import numpy as np
import pickle
import pandas as pd
import time
import os
import init_common as x

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

def fun_log_list(x,colname,fun=np.mean,ma=False,w=6,s=1):
    z = pd.DataFrame(x[0][colname])
    if len(z)>1:
        for idx in x[1:]:
            u = pd.DataFrame(np.log(idx[colname]*s))
            z = pd.concat([z,u],1)
    z_res = z.apply(fun,1)
    if ma:
        z_ma = z_res.rolling(window=w).mean()
        z_ma[0:w-1] = z_res.rolling(window=2).mean()[0:w-1]
        z_res[1:]= z_ma[1:]
    return z_res

def mean_list(x,colname,ma=False,w=6):
    z = pd.DataFrame(x[0][colname])
    if len(z)>1:
        for idx in x[1:]:
            z = pd.concat([z,pd.DataFrame(idx[colname])],1)
    z_mean = z.apply(np.mean,1)
    if ma:
        z_mean0 = z_mean.rolling(window=w).mean()
        z_mean0[0:w-1] = z_mean.rolling(window=2).mean()[0:w-1]
        z_mean[1:]=z_mean0[1:]
    return z_mean

def quant_list(x,colname,val,ma=False,w=6):
    z = pd.DataFrame(x[0][colname])
    if len(z)>1:
        for idx in x[1:]:
            z = pd.concat([z,pd.DataFrame(idx[colname])],1)
    z_quant = z.apply(np.quantile,axis=1,args=(val,))
    if ma:
        z_quant0 = z_quant.rolling(window=w).mean()
        z_quant0[0:w-1] = z_quant.rolling(window=2).mean()[0:w-1]
        z_quant[1:]=z_quant0[1:]
    return z_quant


def std_list(x,colname,ma=False,w=6):
    z = pd.DataFrame(x[0][colname])
    if len(z)>1:
        for idx in x[1:]:
            z = pd.concat([z,pd.DataFrame(idx[colname])],1)  
    z_std = z.apply(np.std,1)        
    if ma:
        z_std0 = z_std.rolling(window=w).mean()
        z_std0[0:w-1] = z_std.rolling(window=2).mean()[0:w-1]
        z_std[1:]=z_std0[1:]
    return z_std

def mean_diff_list(x,colname,w=8):
    u = pd.array(x[0][colname])
    z = pd.DataFrame(u[1:]-u[:-1])
    z_mean = z
    n = len(x)
    if n>1:
        for idx in np.arange(n):
            u = pd.array(x[idx][colname])
            z = pd.concat([z,pd.DataFrame(u[1:]-u[:-1])],1)
        z_mean = z.apply(np.mean,1)
        z_mean0 = z_mean.rolling(window=w).mean()
        z_mean0[0:w-1] = z_mean.rolling(window=2).mean()[0:w-1]
        z_mean[1:]=z_mean0[1:]
    return z_mean


def std_diff_list(x,colname):
    u = pd.array(x[0][colname])
    z = pd.DataFrame(u[1:]-u[:-1])
    n = len(x)
    if n>1:
        for idx in np.arange(n):
            u = pd.array(x[idx][colname])
            z = pd.concat([z,pd.DataFrame(u[1:]-u[:-1])],1)
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

def create_alpha_matrix(idx_sel,init_sel,tar_sel,per,init_rem,tar_rem):
    alpha_mat = np.full([x.N_locs,x.N_per],init_rem)
    x0 = np.full([x.N_locs,1],init_rem)
    x1 = np.full([x.N_locs,1],tar_rem)
    alpha_mat = np.linspace(x0,x1,x.N_per)
    alpha_mat[idx_sel,0:per] = np.linspace(init_sel,tar_sel,per)
    return alpha_mat   
    
def get_paths(R0,alpha):
    out_filename_root = "./out"
    out_fig_root = "./fig"
    out_stat_root = "./stat"
    try:
        R0_str = str(int(100*R0))
    except:
        R0_str = R0
    try:
        alpha_str = str(int(100*alpha))
    except:
        alpha_str = alpha
    foldname = "R0_"+R0_str+"_alpha_"+alpha_str
    out_filename_dir = out_filename_root+"/"+foldname
    out_fig_dir = out_fig_root+"/"+foldname
    out_stat_dir = out_stat_root+"/"+foldname
    return out_filename_dir,out_fig_dir,out_stat_dir 
    
    
def setup_paths(R0,alpha):
    
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
    try:
        R0_str = str(int(100*R0))
    except:
        R0_str = R0
    try:
        alpha_str = str(int(100*alpha))
    except:
        alpha_str = alpha
    foldname = "R0_"+R0_str+"_alpha_"+alpha_str
    out_filename_dir = out_filename_root+"/"+foldname
    out_fig_dir = out_fig_root+"/"+foldname
    out_stat_dir = out_stat_root+"/"+foldname
    try:
        os.mkdir(out_filename_dir)
    except:
        ethrown=True
    try:
        os.mkdir(out_fig_dir)
    except:
        ethrown=True
    try:
        os.mkdir(out_stat_dir)
    except:
        ethrown=True
    return out_filename_dir,out_fig_dir,out_stat_dir 
    
    

