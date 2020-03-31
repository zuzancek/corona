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

def set_params(R0_type):
    if R0_type == 1:
        R0 = 2.2
    elif R0_type == 2:
        R0 = 2.6
    beta = R0*x.gamma
    return beta
