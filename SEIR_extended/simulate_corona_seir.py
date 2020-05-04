import sys
import numpy as np
import pickle
import pandas as pd
import time
import init as x
from random import sample

def dsim():
    # epidemiology
    con_share = 1/3
    first_con = np.round(x.first_infections*con_share)
    first_exp = np.round(x.first_infections*con_share)
    first_rem = x.first_infections-first_con-first_exp
    SECR = np.zeros(shape=(x.N_locs, 4)) 
    SECR[:,0] = x.N_popul- x.first_infections
    SECR[:,1] = first_con
    SECR[:,2] = first_exp
    SECR[:,3] = first_rem
    row_sums = SCR.sum(axis=1)
    SCR_norm = SCR / row_sums[:, np.newaxis]
      
    rho_eff = x.rho_eff 
    S = np.zeros((x.N_per,))
    E = np.zeros((x.N_per,))
    C = np.zeros((x.N_per,))
    R = np.zeros((x.N_per,))
    R0 = x.R0_mean# = np.array(sample(x.R0_vec_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    T0 = x.Tinf_mean
    T1 = x.Tinc_mean
    print([T0,T1,R0])
    
    for t in range(x.N_per):
        c_mat = np.array([SCR_norm[:,1],]*x.N_locs).transpose()
        OD_con = np.round(x.OD*c_mat)#*rho_eff
        c_in = np.round(OD_con.sum(axis=0))#*alpha_vec)
        s_out = R0/T0*SCR[:, 0]*c_in/(x.N_popul + x.OD.sum(axis=0))
        c_out = SCR[:,1]/T1
        SCR[:,0] = SCR[:,0]-s_out
        SCR[:,1] = SCR[:,1]+s_out-c_out  
        SCR[:,2] = SCR[:,2]+c_out          
        row_sums = SCR.sum(axis=1)                
        SCR_norm = SCR/row_sums[:, np.newaxis]
        agg_epi = SCR.sum(axis=0)/x.N_popul_size  
        S[t] = agg_epi[0]
        C[t] = agg_epi[1]
        R[t] = agg_epi[2] #1-agg_epi[0]-agg_epi[1]-agg_epi[2]
        print([S[t],C[t],R[t]])
    res = pd.DataFrame(list(zip(S,C,R)), columns = ['S','C','R'])
    
    return res 
          
        
    
def simul(beta_list,Tinf_inv_list,Tinc_inv_list,alpha_mat,it): 
    
    # epidemiology
    con_share = 1/3
    first_con = np.round(x.first_infections*con_share)
    first_rem = x.first_infections-first_con
    SECR = np.zeros(shape=(x.N_locs, 4)) 
    SECR[:,0] = x.N_popul- x.first_infections
    SECR[:,2] = first_con
    SECR[:,3] = first_rem
    row_sums = SECR.sum(axis=1)
    SECR_norm = SECR / row_sums[:, np.newaxis]
      
    rho_eff = x.rho_eff 
    S = np.zeros((x.N_per,))
    E = np.zeros((x.N_per,))
    C = np.zeros((x.N_per,))
    R = np.zeros((x.N_per,))
    
    beta_mat = np.array(sample(beta_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    Tinc_inv_mat = np.array(sample(x.Tinc_inv_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    Tinf_inv_mat = np.array(sample(x.Tinf_inv_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    rho_eff = x.rho_eff
    Tinc = x.Tinc_mean
    Tinf = x.Tinf_mean
    print([Tinf,Tinc])
    R0_mat = np.array(sample(x.R0_vec_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    
    cont = True
    while cont:
        try:
            for t in range(x.N_per):
                beta_vec = beta_mat[:,t]
                alpha_vec = alpha_mat[:,t]
                Tinc_inv_vec = Tinc_inv_mat[:,t]
                Tinf_inv_vec = Tinf_inv_mat[:,t]
                # Trec_unobs_asymp_vec = Trec_unobs_asymp[:,t]
                # Trec_obs_asymp_vec = Trec_obs_asymp[:,t]
                # matrix of contageous
                con_mat = np.array([SECR_norm[:,2],]*x.N_locs).transpose()
                OD_con = np.round(x.OD*con_mat)#*rho_eff
                ## Epidemiology: use stochastic approach
                # 0: suspectible
                con_in = np.round(OD_con.sum(axis=0))#*alpha_vec)
                exp_in = beta_vec*SECR[:, 0]*con_in/(x.N_popul + x.OD.sum(axis=0))
                exp_in = np.where(exp_in>SECR[:,0],SECR[:,0],exp_in)
                exp_out = Tinc_inv_vec*SECR[:,1]
                con_out = Tinf_inv_vec*SECR[:,2]
                SECR[:,0] = SECR[:,0]-exp_in
                # 1: exposed
                SECR[:,1] = SECR[:,1]+exp_in-exp_out  
                # 2: contageous (infectious)
                SECR[:,2] = SECR[:,2]+exp_out-con_out
                # 3: removed, immune & non contageous  (but still infected)
                SECR[:,3] = SECR[:,3]+con_out            
                ## normalisation, copying
                SECR = np.where(SECR<0,0,SECR)     
                row_sums = SECR.sum(axis=1)                
                SECR_norm = SECR/row_sums[:, np.newaxis]
                agg_epi = SECR.sum(axis=0)/x.N_popul_size  
                S[t] = agg_epi[0]
                E[t] = agg_epi[1]
                C[t] = agg_epi[2]
                R[t] = agg_epi[3] #1-agg_epi[0]-agg_epi[1]-agg_epi[2]
                print([S[t],E[t],C[t],R[t]])
            cont = False
        except: 
            print("exc")
    res = pd.DataFrame(list(zip(S,E,C,R)), columns = ['S','E','C','R'])
    
    return res 
