import sys
import numpy as np
import pickle
import pandas as pd
import time
import init_common as x
from random import sample


def simul(beta_list,Trec_list,alpha_mat): 
    SIR = np.zeros(shape=(x.N_locs, 3)) 
    SIR[:,0] = x.N_popul- x.first_infections
    SIR[:, 1] = SIR[:, 1] + x.first_infections    
    row_sums = SIR.sum(axis=1)
    SIR_n = SIR / row_sums[:, np.newaxis]
    SIR_sim = SIR.copy()
    SIR_nsim = SIR_n.copy()
    
    public_trans_vec = np.full(x.N_locs, public_trans) ## <--- change this
    infected_pop_norm = np.zeros((x.N_per,))
    susceptible_pop_norm = np.zeros((x.N_per,))
    recovered_pop_norm = np.zeros((x.N_per,))
    SIR_sim_arr=np.zeros((SIR_sim.shape[0],SIR_sim.shape[1],x.N_per)) 
    beta_mat = np.array(sample(beta_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    rec_mat = np.array(sample(Trec_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    gamma_mat = 1.0/rec_mat
    cont = True
    while cont:
        try:
            for j in (range(x.N_per)): 
                beta_vec = beta_mat[:,j]
                gamma_vec = gamma_mat[:,j]
                rec_time_vec = rec_mat[:,j]
                alpha_vec = alpha_mat[:,j]
                # Matice infekcii
                infected_mat = np.array([SIR_nsim[:,1],]*x.N_locs).transpose()
                OD_infected = np.round(x.OD*infected_mat)
                inflow_infected = OD_infected.sum(axis=0)
                inflow_infected = np.round(inflow_infected*alpha_vec)
                new_infect = beta_vec*SIR_sim[:, 0]*inflow_infected/(x.N_popul + x.OD.sum(axis=0))
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
        except: 
            print("exc")
        
    res = pd.DataFrame(list(zip(infected_pop_norm, susceptible_pop_norm, recovered_pop_norm)), columns = ['inf','sus','rec'])
    
    return res 
