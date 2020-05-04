import sys
import numpy as np
import pickle
import pandas as pd
import time
import init as x
from random import sample
import init_healthcare as xx


def get_htc(inflow): 
    
    HTC = np.zeros(shape=(x.N_per+1, 13)) 
    HTC[0,0] = x.first_infections_unobserved_asymptomatic
    HTC[0,1] = x.first_infections_unobserved_symptomatic
    HTC[0,2] = x.first_infections_observed_asymptomatic
    HTC[0,3] = x.first_infections_mild
    HTC[0,4] = x.first_infections_hospital
    HTC[0,5] = x.first_infections_icu 
    
    omega_obs = xx.omega_obs
    omega_unobs_asymp = xx.omega_asymp
    omega_obs_asymp = xx.omega_asymp
    prob_hosp_obs = xx.prob_hosp_obs/(1-omega_obs_asymp)
    prob_hosp_unobs = xx.prob_hosp_unobs/(1-omega_unobs_asymp)
    prob_icu = xx.prob_icu
    prob_vent = xx.prob_vent
    prob_surv = xx.prob_surv  
    T_death = xx.T_death-xx.T_hosp
    Trec_unobs_asymp = xx.Trec_asymp
    Trec_obs_asymp = xx.Trec_asymp
    Trec_unobs_symp = xx.Trec_mild
    Trec_mild = xx.Trec_mild
    Trec_hosp = xx.Trec_hosp-xx.T_hosp
    Trec_icu = xx.Trec_icu-xx.T_hosp
    
    for t in range(x.N_per):
        # 0: unobserved, asymptomatic
        inf_unobs_asymp_in = (1-omega_obs)*omega_unobs_asymp*inflow[t]
        inf_unobs_asymp_out = HTC[t,0]/Trec_unobs_asymp
        HTC[t+1,0] = HTC[t,0]+inf_unobs_asymp_in-inf_unobs_asymp_out
        # 1: unobserved, symptomatic (in fact mild cases only)               
        inf_unobs_symp_in = (1-omega_obs)*(1-omega_unobs_asymp)*inflow[t]
        inf_unobs_symp_out = (prob_hosp_unobs+1/Trec_unobs_symp)*HTC[t,1]
        HTC[t+1,1] = HTC[t,1]+inf_unobs_symp_in-inf_unobs_symp_out               
        # 2: observed, asymptomatic
        inf_obs_asymp_in = (omega_obs)*omega_obs_asymp*inflow[t]
        inf_obs_asymp_out = HTC[t,2]/Trec_obs_asymp
        HTC[t+1,2] = HTC[t,2]+inf_obs_asymp_in-inf_obs_asymp_out
        # 3: observed, mild cases
        inf_obs_mild_in = omega_obs*(1-omega_obs_asymp)*inflow[t]
        inf_obs_mild_out = (prob_hosp_obs+1/Trec_mild)*HTC[t,3]
        HTC[t+1,3] = HTC[t,3]+inf_obs_mild_in-inf_obs_mild_out 
        # 4: serious/hospital cases
        inf_hosp_in = prob_hosp_obs*HTC[t,3]+prob_hosp_unobs*HTC[t,1]
        inf_hosp_out = (prob_icu+1/Trec_hosp)*HTC[t,4]
        HTC[t+1,4] = HTC[t,4]+inf_hosp_in-inf_hosp_out
        # 5: ICU cases
        inf_icu_in = prob_icu*HTC[t,4]
        inf_icu_out = (prob_surv+1/T_death)*HTC[t,5]
        HTC[t+1,5] = HTC[t,5]+inf_icu_in-inf_icu_out                
        # 6: Recovered, unobserved
        rec_unobs_in = inf_unobs_asymp_out+1/Trec_unobs_symp*HTC[t,1]
        HTC[t+1,6] = HTC[t,6]+rec_unobs_in
        # 7: Recovered, observed
        rec_obs_in = inf_obs_asymp_out+1/Trec_mild*HTC[t,3]+1/Trec_hosp*HTC[t,4]+prob_surv*HTC[t,5]
        HTC[t+1,7] = HTC[t,7]+rec_obs_in
        # 8: Dead (observed only)
        HTC[t+1,8] = HTC[t,8]+1/T_death*HTC[t,5]  
        # 9: Observed
        HTC[t+1,9] = HTC[t+1,2]+HTC[t+1,3]+HTC[t+1,4]+HTC[t+1,5]
        # 10: UnObserved
        HTC[t+1,10] = HTC[t+1,0]+HTC[t+1,1]
        # 11: Ventilation
        HTC[t+1,11] = prob_vent*HTC[t+1,5]
        # 12: Hospital
        HTC[t+1,12] = HTC[t+1,5]+HTC[t+1,4]
    
    return HTC 

