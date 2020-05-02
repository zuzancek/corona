import sys
import numpy as np
import time
import os

def get_healthcare_params(type):
    if type == 0: # IZP
        repinf_inf_ratio = 1
        symp_repinf_ratio = .6
        hosp_symp_ratio = .074306
        icu_hosp_ratio = .1375
    return repinf_inf_ratio,symp_repinf_ratio,hosp_symp_ratio,icu_hosp_ratio
