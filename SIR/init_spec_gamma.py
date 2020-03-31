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
gamma_scale = 2
rec_vec = np.random.gamma(Trec/gamma_scale,gamma_scale,10000000)
Trec_mean = np.average(rec_vec)
gamma_mean = 1/Trec_mean  
gamma = 0.1

beta = x.R0*gamma
beta_scale = 1/10

