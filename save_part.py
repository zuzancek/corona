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
import init as x
import os

def save_file_pickle(data,filenameroot):
    
    fid = open(filenameroot+'.pickle','wb')
    N = len(data)
    idx = N-1
    idx0 = idx
    idx1 = 0
    cont = True
    k = 0
    while cont:
        y = data[0:idx]
        try:
            pickle.dump(y,fid)
            excthrown = False
        except:
            excthrown = True
        if excthrown:
            idx0 = idx
            if idx<=1: 
                cont = False
            else: 
                if idx1>0: 
                    idx = int((idx+idx1)/2)
                else: 
                    idx = int((idx+1)/2)-1
        else:
            if idx1==idx or idx>=N-2: 
                cont = False
            else: 
                idx = int((idx+idx0)/2)
                idx1 = idx
        
    try:
        fid.close()
    except:
        print('file closed')
        
    filepartsize = idx+1
    filecnt = 1+int(N/filepartsize)
    if filecnt>1:
        k=0
        for fidx in range(filecnt):
            filename0 = filenameroot+'_'+str(k)+'.pickle'
            fid0 = open(filename0,'wb')
            lastidx = min((k+1)*filepartsize,N)
            y = data[k*filepartsize:lastidx]
            pickle.dump(y,fid0)
            k+=1
            fid0.close()
