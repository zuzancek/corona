import sys
import numpy as np
import pickle
import pandas as pd
import time
import init_common as x
import init_spec_inout as xx
from random import sample


def simul(beta0_list,beta1_list,Tinf_list,dump_vect0,dump_vect1,it): 
    
    SIR = np.zeros(shape=(x.N_locs, 3)) 
    SIR[:,0] = x.N_popul- x.first_infections
    SIR[:,1] = SIR[:, 1] + x.first_infections    
    row_sums = SIR.sum(axis=1)
    SIR_n = SIR / row_sums[:, np.newaxis]
    SIR_sim = SIR.copy()
    SIR_nsim = SIR_n.copy()
    np.fill_diagonal(x.OD,0)
    OD = xx.OD
    
    rho_eff = x.rho_eff
    I_norm = np.zeros((x.N_per,))
    S_norm = np.zeros((x.N_per,))
    R_norm = np.zeros((x.N_per,))
    SIR_sim_arr=np.zeros((SIR_sim.shape[0],SIR_sim.shape[1],x.N_per)) 
    beta0_mat = np.array(sample(beta0_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    beta1_mat = np.array(sample(beta1_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    gamma_mat = 1.0/np.array(sample(Tinf_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    cont = True
    Inew_v = np.zeros((x.N_per,))
    # OD matrix (i,j): from j to i
    # OD : sum(axis=0): all "to" (sum of all row data in a column), sum(axis=1): all "from"
    while cont:
        if True:
            for j in (range(x.N_per)): 
                dump0 = dump_vect0[j]
                dump1 = dump_vect1[j]   
                beta0_vec = dump0*beta0_mat[:,j]
                beta1_vec = dump1*beta1_mat[:,j]
                gamma_vec = gamma_mat[:,j]
                # relative shares of S,I
                y = SIR_sim[:,0]/x.N_popul
                z = SIR_sim[:,1]/x.N_popul
                # term1: infected within their municipality by other inhabitants of the municipality, not during working hours
                out_work = rho_eff*beta0_vec*y*SIR_sim[:,1]                
                # term2: infected within their municipality by comuters from differnt munic., during working hours
                denom = x.N_popul+OD.sum(axis=1)-OD.sum(axis=0)
                # in_work_1 = np.zeros(x.N_locs)
                in_work_1 = (SIR_sim[:,0]-y*OD.sum(0))*(np.sum(OD*(z*beta1_vec),1)+(SIR_sim[:,1]-z*OD.sum(0))*beta0_vec)
                in_work_1 = in_work_1/denom
                # term3               
                # in_work_2 = np.zeros(x.N_locs)
                in_work_2_nom = (SIR_sim[:,1]-z*OD.sum(0))*beta1_vec+np.sum(OD*beta1_vec*z,1)
                in_work_2 = y*np.sum(x.OD.transpose()*in_work_2_nom/denom,1)
                I_new = x.tau*out_work+(1-x.tau)*(in_work_1+in_work_2)
                # S cannot be negative...
                I_new = np.where(I_new>SIR_sim[:,0],SIR_sim[:,0],I_new)
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
                S = SIR_sim[:,0].sum()/x.N_popul_size
                I = SIR_sim[:,1].sum()/x.N_popul_size
                R = SIR_sim[:,2].sum()/x.N_popul_size
                I_norm[j]=I
                S_norm[j]=S
                R_norm[j]=R
                Inew_v[j] = I_new.sum()
            cont = False
        #except: 
            print("exc")
        
    res = pd.DataFrame(list(zip(S_norm,I_norm,R_norm,Inew_v)), columns = ['sus','inf','rec','i_new'])
    
    return res 

    return res 
