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
import init as x

def set_params():
    if x.R0_type == 1:
        x.R0 = 2.0
    elif x.R0_type == 2:
        x.R0 = 2.6
    x.beta = x.R0*x.gamma
