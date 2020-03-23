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
import os, pathlib
import fnmatch as fn
import init as x

def load_file_pickle(filepath,pattern):
    
    # get all files in a given directory
    path = pathlib.Path(filepath)    
    allfiles = []
    reqfiles = []
    print(path)
    filepattern = str(path)+"*"
    print(filepattern)
    for fitem in path.iterdir():
        allfiles.append(fitem)
        if fn.fnmatch(fitem,filepattern):
            reqfiles.append(fitem)
    
    # load files
    filecnt = len(reqfiles)
    if filecnt==0:
        data = []
    elif filecnt==1:
        fid = open(reqfiles[0],'rb')
        data = pickle.load(fid)
        fid.close()        
    else:
        fid = open(reqfiles[0],'rb')
        y = pickle.load(fid)
        data = np.full(len(y)*filecnt, 0)
        fid.close()
        k = 0
        for fitem in reqfiles:
            fid = open(fitem,'rb')
            y = pickle.load(fid)
            data[idx:idx+len(y)] = y
            fid.close()
            idx=+len(y)
        data = data[0:idx]
        
    return data

        