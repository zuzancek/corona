import sys
import numpy as np
import time
import os
import init_common as x

def get_healthcare_params(type):
    if type == 0: # IZP
        repinf_inf_ratio = 1
        symp_repinf_ratio = .6
        hosp_symp_ratio = .074306
        icu_hosp_ratio = .1375
    elif type == 1: # sth.else
        repinf_inf_ratio = 1/x.first_infections_correction_multiplier
        symp_repinf_ratio = 1
        hosp_symp_ratio = .2/symp_repinf_ratio
        icu_hosp_ratio = .2
        
    return repinf_inf_ratio,symp_repinf_ratio,hosp_symp_ratio,icu_hosp_ratio
