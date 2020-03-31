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
beta = x.R0*x.gamma
beta_scale = 1/25
gamma_scale = 2

