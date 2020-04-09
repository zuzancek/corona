import sys
import numpy as np
import pickle
import pandas as pd
import time
import init_common as x
from random import sample


def simul(beta_list,Trec_list,alpha_mat,it): 
    
    SIR = np.zeros(shape=(x.N_locs, 9)) 
    # S
    # I_u, I_a, I_m, I_s,I_c,
    # R_u,R_o,D
    SIR[:,0] = x.N_popul- x.first_infections
    SIR[:,1] = SIR[:,1]+x.first_infections_unobserved
    SIR[:,2] = SIR[:,2]+x.first_infections_asymptomatic
    SIR[:,3] = SIR[:,3]+x.first_infections_mild
    SIR[:,4] = SIR[:,4]+x.first_infections_serious
    SIR[:,5] = SIR[:,5]+x.first_infections_icu 
    row_sums = SIR.sum(axis=1)
    SIR_norm = SIR / row_sums[:, np.newaxis]
    #SIR_sim = SIR.copy()
    #SIR_norm = SIR_n.copy()
    idx_i_0 = 1
    idx_i_1 = 5
    
    rho_trans_vec = x.rho_trans_vec #5x1
    S = np.zeros((x.N_per,))
    I_u = np.zeros((x.N_per,))
    I_a = np.zeros((x.N_per,))
    I_m = np.zeros((x.N_per,))
    I_s = np.zeros((x.N_per,))
    I_c = np.zeros((x.N_per,))
    I_o = np.zeros((x.N_per,))
    R_o = np.zeros((x.N_per,))
    R_u = np.zeros((x.N_per,))
    R = np.zeros((x.N_per,))
    I = np.zeros((x.N_per,))
    D = np.zeros((x.N_per,))
    
    #SIR_sim_arr=np.zeros((SIR_sim.shape[0],SIR_sim.shape[1],x.N_per)) 
    beta_mat = np.array(sample(beta_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    rec_mat = np.array(sample(Trec_list,x.N_per*x.N_locs)).reshape((x.N_locs,x.N_per))
    gamma_mat = 1.0/rec_mat
    cont = True
    while cont:
        try:
            for j in (range(x.N_per)): 
                beta_vec = beta_mat[:,j]
                alpha_vec = alpha_mat[:,j]
                # Matice infekcii
                inf_mat = np.array([SIR_norm[:,idx_i_0:idx_i_1],]*x.N_locs).transpose()
                OD_inf = np.round(x.OD*inf_mat)*rho_trans_vec
                inf_in = OD_inf.sum(axis=0)
                inf_in = np.round(inf_in*alpha_vec)
                inf_new = beta_vec*SIR[:, 0]*inflow_infected/(x.N_popul + x.OD.sum(axis=0))
                SIR[:,0] = SIR[:,0]-inf_new
                inf_unobs_in = delta_unobs*inf_new
                inf_unobs_out = (phi_unobs+psi_unobs)*SIR[:,1]
                SIR[:,1] = SIR[:,1]+inf_unobs_in-inf_unobs_out
                inf_asymp_in = (1-delta_unobs)*delta_asymp*inf_new
                inf_asymp_out = (phi_asymp+psi_asymp)*SIR[:,2]
                SIR[:,2] = SIR[:,2]+inf_asymp_in-inf_asymp_out
                inf_mild_in = (1-delta_unobs)*(1-delta_asymp)*inf_new+psi_asymp*SIR[:,2]
                inf_mild_out = (phi_mild+psi_mild)*SIR[:,3]
                SIR[:,3] = SIR[:,3]+inf_mild_in-inf_mild_out
                inf_sev_in = psi_mild*SIR[:,3]+psi_unobs*SIR[:,1]
                inf_sev_out = (phi_sev+psi_sev)*SIR[:,4]
                SIR[:,4] = SIR[:,4]+inf_sev_in-inf_sev_out
                inf_icu_in = psi_sev*SIR[:,4]
                inf_icu_out = (phi_icu+psi_sev)*SIR[:,4]
                SIR[:,5] = SIR[:,5]+inf_icu_in-inf_icu_out
                SIR[:,6] = SIR[:,6]+phi_unobs*SIR[:,1]
                SIR[:,7] = SIR[:,7]+phi_asymp*SIR[:,2]+phi_mild*SIR[:3]+phi_sev*SIR[:,4]+phi_icu*SIR[:,5]
                SIR[:,8] = SIR[:,8]+psi_sev*SIR[:,4]
                SIR = np.where(SIR<0,0,SIR)
                row_sums = SIR.sum(axis=1)                
                SIR_norm = SIR/row_sums[:, np.newaxis]
                agg = SIR.sum(axis=0)/x.N_popul_size
                S[j] = agg[0]
                I_u[j] = agg[1]
                I_a[j] = agg[2]
                I_m[j] = agg[3]
                I_s[j] = agg[4]
                I_c[j] = agg[5]
                I_o[j] = I_a[j]+I_m[j]6I_s[j]+I_c[j]
                I[j] = I_u[j]+I_o[j]
                R_u[j] = agg[6]
                R_o[j] = agg[7]
                D[j] = agg[8]
                R[j] = R_u[j]+R_o[j]+D[j]
            cont = False
        except: 
            print("exc")
        
    res = pd.DataFrame(list(zip(S,I,I_u,I_o,I_a,I_m,I_s,I_c,R_u,R_o,D,R)), columns = ['sus','inf','inf_unobs','inf_obs','inf_asymp','inf_mild','inf_sev','inf_icu','rec_unobs','rec_obs','death','rec'])
    
    return res 
