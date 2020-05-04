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
import os
import helpers as hp
from random import sample
from sklearn.utils import shuffle
import init

## clinical parameters:
# probabilities & shares
omega_obs = 1.0/init.first_infections_correction_multiplier
omega_asymp = .5
prob_icu = 1/6
prob_surv = 0.95
prob_vent = 0.4
# times
T_hosp = 5
Trec_hosp = T_hosp+8
Trec_icu = T_hosp+10
T_death = 15+T_hosp
Trec_mild = 8
Trec_asymp = 4
prob_hosp = .2*Trec_asymp/(Trec_hosp)*(1-omega_asymp)*omega_obs
prob_hosp_obs = .2*Trec_asymp/(Trec_hosp)*(1-omega_asymp)
prob_hosp_unobs = .02*Trec_asymp/(Trec_hosp)*(1-omega_asymp)*omega_obs

first_infections_original_num = sum(init.first_infections_original)
first_infections_observed_asymptomatic = int(init.first_infections_original_num*omega_asymp)
first_infections_observed_symptomatic = init.first_infections_original_num-first_infections_observed_asymptomatic
first_infections_mild = first_infections_observed_symptomatic
first_infections_hospital = int(1/3*first_infections_observed_symptomatic)
first_infections_icu = 1
first_infections_unobserved = sum(init.first_infections)-init.first_infections_original_num
first_infections_unobserved_asymptomatic = int(first_infections_unobserved*omega_asymp)
first_infections_unobserved_symptomatic = first_infections_unobserved-first_infections_unobserved_asymptomatic