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
import init_common as x
import init_spec_gamma as xx

def simul(public_trans):
    
    SIR = np.zeros(shape=(x.N_locs_nosen, 3)) 
    SIR[:,0] = x.N_locs_nosen - x.first_infections
    SIR[:, 1] = SIR[:, 1] + x.first_infections    
    
    row_sums = SIR.sum(axis=1)
    SIR_n = SIR / row_sums[:, np.newaxis]
    gamma_vec = np.full(x.N_locs_nosen, x.gamma)
    public_trans_vec = np.full(x.N_locs_nosen, public_trans)
    SIR_sim = SIR.copy()
    SIR_nsim = SIR_n.copy()
    infected_pop_norm = np.zeros((x.N_simul,))
    susceptible_pop_norm = np.zeros((x.N_simul,))
    recovered_pop_norm = np.zeros((x.N_simul,))
    SIR_sim_arr=np.zeros((SIR_sim.shape[0],SIR_sim.shape[1],x.N_simul)) 
    beta_mat = np.random.gamma(xx.beta/xx.beta_scale,xx.beta_scale, size=[x.N_locs_nosen,x.N_simul])
    cont = True
    while cont:
        try:
            for j in (range(x.N_per)):        
                beta_vec = beta_mat[:,j]
                # Matice infekcii
                infected_mat = np.array([SIR_nsim[:,1],]*x.N_locs_nosen).transpose()
                OD_infected = np.round(x.OD*infected_mat)
                inflow_infected = OD_infected.sum(axis=0)
                inflow_infected = np.round(inflow_infected*public_trans_vec)
                new_infect = beta_vec*SIR_sim[:, 0]*inflow_infected/(x.N_popul_nosen + x.OD.sum(axis=0))
                new_recovered = gamma_vec*SIR_sim[:, 1]
                new_infect = np.where(new_infect>SIR_sim[:, 0], SIR_sim[:, 0], new_infect)
                SIR_sim[:, 0] = SIR_sim[:, 0] - new_infect
                SIR_sim[:, 1] = SIR_sim[:, 1] + new_infect - new_recovered
                SIR_sim[:, 2] = SIR_sim[:, 2] + new_recovered
                SIR_sim = np.where(SIR_sim<0,0,SIR_sim)
                row_sums = SIR_sim.sum(axis=1)
                SIR_nsim = SIR_sim / row_sums[:, np.newaxis]
                SIR_sim_arr[:,:,j]=SIR_sim
                S = SIR_sim[:,0].sum()/x.N_popul_size
                I = SIR_sim[:,1].sum()/x.N_popul_size
                R = SIR_sim[:,2].sum()/x.N_popul_size
                infected_pop_norm[j] = I
                susceptible_pop_norm[j] = S
                recovered_pop_norm[j] = R
            cont = False
            # print("suc")
        except: 
            print("exc")
        
    ## Vytvor konecnu maticu
    res = pd.DataFrame(list(zip(infected_pop_norm, susceptible_pop_norm, recovered_pop_norm)), columns = ['inf','sus','rec'])
    
    return res,SIR_sim_arr
