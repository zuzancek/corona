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
    I_norm = np.zeros((x.N_per,))
    S_norm = np.zeros((x.N_per,))
    R_norm = np.zeros((x.N_per,))
    SIR_sim_arr=np.zeros((SIR_sim.shape[0],SIR_sim.shape[1],x.N_per)) 
    beta_mat = np.array(sample(beta_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    rec_mat = np.array(sample(Trec_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    gamma_mat = 1.0/rec_mat
    cont = True
    # OD matrix (i,j): from j to i
    # OD : sum(axis=0): all "to" (sum of all row data in a column), sum(axis=1): all "from"
    while cont:
        try:
            for j in (range(x.N_per)): 
                beta_vec = beta_mat[:,j]
                gamma_vec = gamma_mat[:,j]
                rec_time_vec = rec_mat[:,j]
                alpha_vec = alpha_mat[:,j]
                # relative shares of S,I
                y = SIR_sim[:,0]/x.N_popul
                x = SIR_sim[:,1]/x.N_popul
                # term1: infected within their municipality by other inhabitants of the municipality, not during working hours
                out_work = beta_vec*y*SIR_sim[:,1]
                # term2: infected within their municipality by comuters from differnt munic., during working hours
                denom = x.N_popul+OD.sum(axis=1)-OD.sum(axis=0)
                in_work_1 = (SIR_sim[:,0]-y*OD.sum(0))*((SIR_sim[:,1]-x*OD.sum(0))*beta_vec+np.sum(OD*beta_vec*x),1))
                # term3
                in_work_2_nom = (SIR_sim[:,1]-x*OD.sum(0))*beta_vec+np.sum(OD*beta_vec*x,1)
                in_work_2 = y*np.sum(OD.transpose()*in_work_2_nom/denom,1)
                I_new = x.tau*out_work+alpha_vec*(1-x.tau)*(in_work_1+in_work_2)
                # S cannot be negative...
                I_new = np.where(new_I>SIR_sim[:,0],SIR_sim[:0],I_new)
                # Recovered
                R_new = gamma_vec*SIR_sim[:,1]
                # apply changes
                SIR_sim[:,0] = SIR_sim[:,0]-I_new
                SIR_sim[:,1] = SIR_sim[:,1]+I_new-R_new
                SIR_sim[:,2] = SIR_sim[:,2]+R_new
                # correction for negative values
                SIR_sim = np.where(SIR_sim<0,0,SIR_sim)
                # normalise
                rsum = SIR_sim.sum(axis=1)
                S = SIR_sim[:,0].sum()/N_popul_size
                I = SIR_sum[:,1].sum()/N_popul_size
                R = SIR_sum[:,2].sum()/N_popul_size
                I_norm.append(I)
                S_norm.append(S)
                R_norm.append(R)
            cont = False
        except: 
            print("exc")
        
    res = pd.DataFrame(list(zip(S_norm,I_norm,R_norm)), columns = ['sus','inf','rec'])
    
    return res 
