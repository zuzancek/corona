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
sns.set(rc={'figure.figsize':(11, 4)})

## Funkcia pre výpočet priemeru zo simulácií
def sumlist(x):
    tmp=x[0]
    for i in x[1:]:
        tmp=tmp+i
    return tmp/len(x)

## Údaje k počtu obyvateľov na obec 
pop = pd.read_excel('./src/munic_pop.xlsx')
pop_N = np.array(pop['popul'])

## Priradenie GPS suradnic pre kazdu obec
def get_coors_long(x):
    return float(df_coords.loc[df_coords.IDN4.apply(str)==x,'long'])

def get_coors_lat(x):
    return float(df_coords.loc[df_coords.IDN4.apply(str)==x,'lat'])
df_coords=pd.read_excel('./src/obce1.xlsx')
data_i=pop
data_i.loc[:,'long']=data_i.munic.apply(str).apply(get_coors_long)
data_i.loc[:,'lat']=data_i.munic.apply(str).apply(get_coors_lat)

with open('./src/OD_final_16.3.2020.pickle','rb') as f:
    OD=pickle.load(f)


nakazy_sk=pd.DataFrame({'kod':[529346,529346,529320 , 512036,508063,500011,506338,598186,508438,506745,508217,503011,
                    501433,527106,505315,505315,517461,505820,526355],
                      'pocet':[4,3,25,7,4,3,3,2,2,2,2,2,1,1,1,1,1,1,1]})

first_infections=np.zeros(2926)
for i in np.arange(nakazy_sk.shape[0]):
    first_infections[pop.munic==nakazy_sk.kod.iloc[i]]=nakazy_sk.pocet.iloc[i]
    
first_infections_original=first_infections
first_infections=first_infections_original*6

## definition of key parameters
beta = 0.24 # "Transmission rate" <------ TEST HERE 2.4 with confint <2.2,2.6>
gamma = 0.10 # "Recovery rate", length of sickness is 12 days approx.
R0 = beta/gamma # Reprodukcne cislo ("Basic reproduction number") 

## technical params
N_k = pop.popul.to_numpy()          # Populacia
locs_len = len(N_k)                 # Pocet obci
simul_len = 200   
simul_cnt = 50
public_trans_glob = 0.5

data_senior=pd.read_excel('./src/senior.xlsx')
data_senior.loc[:,'munic']=data_senior.munic.apply(lambda x: x[-6:]).apply(int)
data_senior=data_senior.sort_values(by=['munic'])
