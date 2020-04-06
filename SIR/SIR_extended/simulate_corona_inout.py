import sys
import numpy as np
import pickle
import pandas as pd
import time
import init_common as x
from random import sample


def simul(beta_list,Trec_list,alpha_in_mat,alpha_out_mat,it): 
    
    SIR = np.zeros(shape=(x.N_locs, 3)) 
    SIR[:,0] = x.N_popul- x.first_infections
    SIR[:,1] = SIR[:, 1] + x.first_infections    
    row_sums = SIR.sum(axis=1)
    SIR_n = SIR / row_sums[:, np.newaxis]
    SIR_sim = SIR.copy()
    SIR_nsim = SIR_n.copy()
    np.fill_diagonal(x.OD,0)
    OD=x.OD
    
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
        if True:
            for j in (range(x.N_per)): 
                beta_vec = beta_mat[:,j]
                gamma_vec = gamma_mat[:,j]
                rec_time_vec = rec_mat[:,j]
                alpha_in_vec = alpha_in_mat[:,j]
                alpha_out_vec = alpha_out_mat[:,j]                
                print(1)
                # relative shares of S,I
                y = SIR_sim[:,0]/x.N_popul
                z = SIR_sim[:,1]/x.N_popul
                # term1: infected within their municipality by other inhabitants of the municipality, not during working hours
                out_work = beta_vec*y*SIR_sim[:,1]                
                # term2: infected within their municipality by comuters from differnt munic., during working hours
                denom = x.N_popul+OD.sum(axis=1)-OD.sum(axis=0)
                print(2)
                in_work_1 = np.zeros*(x.N_locs)
                in_work_1 = (SIR_sim[:,0]-y*OD.sum(0))*((SIR_sim[:,1]-z*OD.sum(0))*beta_vec+np.sum(OD*(z*beta_vec),1))
                in_work_1 = in_work_1/denom
                # term3
                
                print(3)
                in_work_2_nom = (SIR_sim[:,1]-z*OD.sum(0))*beta_vec+np.sum(OD*beta_vec*z,1)
                in_work_2 = y*np.sum(x.OD.transpose()*in_work_2_nom/denom,1)
                I_new = x.tau*alpha_in_vec*out_work+alpha_out_vec*(1-x.tau)*(in_work_1+in_work_2)
                # S cannot be negative...
                I_new = np.where(new_I>SIR_sim[:,0],SIR_sim[:0],I_new)
                # Recovered
                R_new = gamma_vec*SIR_sim[:,1]
                
                print(4)
                # apply changes
                SIR_sim[:,0] = SIR_sim[:,0]-I_new
                SIR_sim[:,1] = SIR_sim[:,1]+I_new-R_new
                SIR_sim[:,2] = SIR_sim[:,2]+R_new
                # correction for negative values
                SIR_sim = np.where(SIR_sim<0,0,SIR_sim)
                # normalise
                rsum = SIR_sim.sum(axis=1)
                S = SIR_sim[:,0].sum()/x.N_popul_size
                I = SIR_sum[:,1].sum()/x.N_popul_size
                R = SIR_sum[:,2].sum()/x.N_popul_size
                I_norm.append(I)
                S_norm.append(S)
                R_norm.append(R)
                print(j)
            cont = False
        #except: 
            print("exc")
        
    res = pd.DataFrame(list(zip(S_norm,I_norm,R_norm)), columns = ['sus','inf','rec'])
    
    return res 
