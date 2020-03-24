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
    reqfiles = []
    
    filepattern = str(path)+"/"+str(pattern)+"*"
    for fitem in path.iterdir():
        if fn.fnmatch(fitem,filepattern):
            reqfiles.append(fitem)
    
    # load files
    filecnt = len(reqfiles)
    print(filecnt)
    if filecnt==0:
        data = []
    elif filecnt==1:
        fid = open(reqfiles[0],'rb')
        data = pickle.load(fid)
        fid.close()        
    else:
        data = []
        for fitem in reqfiles:
            fid = open(fitem,'rb')
            y = pickle.load(fid)
            data = data+y            
            fid.close()
        
    return data

        