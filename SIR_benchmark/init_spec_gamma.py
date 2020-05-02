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
import init_common as x

# global beta
Trec = x.Tinc+x.Tinf
Trec_mean  = 6.5 # serial interval
Trec_sig = 0.62
m0 = Trec_mean*(Trec_sig**2)
s0 = 1/(Trec_sig**2)
#gamma_scale = 1.25
#rec_vec = np.random.gamma(gamma_scale,Trec/gamma_scale,10000000)
rec_vec = np.random.gamma(m0,s0,10000000)
Trec_mean = np.average(rec_vec)
gamma_mean = 1/Trec_mean  
gamma = 0.1
gamma_vec = 1/rec_vec
gamma_list = gamma_vec.tolist()

beta_mean = x.R0*gamma_mean
beta_scale = 1/10

