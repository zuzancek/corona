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
import helpers as hp

N_nodes = 1000000

## self-isolated
isol_share = 0.7
N_nodes_isol = round(isol_share*N_nodes)
# their beta follows gamma dist. with small variance
beta_scale = 30
beta_isol_mean = 0.15
beta_isol = beta_scale*beta_isol_mean
beta_vec_isol = np.random(beta_isol,1.0/beta_scale,N_nodes_isol)

## non-isolated
a0 = 0.165
a1 = 20
scale_nonisol = 2.75
N_nodes_nonisol = N_nodes-N_nodes_isol
beta_vec_nonisol = hp.power_law_pdf(np.random.uniform(0,1,N_nodes_nonisol),a0,a1,scale_nonisol)

## joint distribution
beta_vec = shuffle(np.transpose([*np.transpose(beta_vec_isol),*np.transpose(beta_vec_nonisol)]),random_state=0)
beta_vec_0 = beta_vec.to_list()
# global beta
beta = x.R0*x.gamma
beta_scale = 1/25
gamma_scale = 2


