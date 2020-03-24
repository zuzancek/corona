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
import init as x

def simul(public_trans):
    
    SIR = np.zeros(shape=(x.locs_len_s, 3)) 
    SIR[:,0] = x.N_k_s                     

    SIR[:, 0] = SIR[:, 0] - x.first_infections
    SIR[:, 1] = SIR[:, 1] + x.first_infections    
    
    row_sums = SIR.sum(axis=1)
    SIR_n = SIR / row_sums[:, np.newaxis]
    gamma_vec = np.full(x.locs_len_s, x.gamma)
    beta_vec = np.full(x.locs_len_s, x.beta)
    public_trans_vec = np.full(x.locs_len_s, public_trans)
    SIR_sim = SIR.copy()
    SIR_nsim = SIR_n.copy()
    infected_pop_norm = np.zeros((x.simul_len,))
    susceptible_pop_norm = np.zeros((x.simul_len,))
    recovered_pop_norm = np.zeros((x.simul_len,))
    SIR_sim_arr=np.zeros((SIR_sim.shape[0],SIR_sim.shape[1],x.simul_len)) 
    cont = True
    while cont:
        j = 0
        try:
            for time_step in (range(x.simul_len)):        
                # beta_vec = beta_mat[:,j] #np.random.gamma(beta, 2, locs_len)
                # Matice infekcii
                infected_mat = np.array([SIR_nsim[:,1],]*x.locs_len_s).transpose()
                OD_infected = np.round(x.OD*infected_mat)
                inflow_infected = OD_infected.sum(axis=0)
                inflow_infected = np.round(inflow_infected*public_trans_vec)
                new_infect = beta_vec*SIR_sim[:, 0]*inflow_infected/(x.N_k_s + x.OD.sum(axis=0))
                new_recovered = gamma_vec*SIR_sim[:, 1]
                new_infect = np.where(new_infect>SIR_sim[:, 0], SIR_sim[:, 0], new_infect)
                SIR_sim[:, 0] = SIR_sim[:, 0] - new_infect
                SIR_sim[:, 1] = SIR_sim[:, 1] + new_infect - new_recovered
                SIR_sim[:, 2] = SIR_sim[:, 2] + new_recovered
                SIR_sim = np.where(SIR_sim<0,0,SIR_sim)
                row_sums = SIR_sim.sum(axis=1)
                SIR_nsim = SIR_sim / row_sums[:, np.newaxis]
                SIR_sim_arr[:,:,j]=SIR_sim
                S = SIR_sim[:,0].sum()/x.N_k.sum()
                I = SIR_sim[:,1].sum()/x.N_k.sum()
                R = SIR_sim[:,2].sum()/x.N_k.sum()
                infected_pop_norm[j] = I
                susceptible_pop_norm[j] = S
                recovered_pop_norm[j] = R
                j=j+1
            cont = False
            # print("suc")
        except: 
            print("exc")
        j+=1
        
    ## Vytvor konecnu maticu
    res = pd.DataFrame(list(zip(infected_pop_norm, susceptible_pop_norm, recovered_pop_norm)), columns = ['inf','sus','rec'])
    
    return res,SIR_sim_arr
