%% initialization & cleanup
initialize;

%% setup
endFcast = dd(2020,12,31); % end of the year
startWave = dd(2020,09,01); % start of the second wave

%% loading stuff
x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility.csv','dateFormat','yyyy-mm-dd','freq','daily');
data = load('results_Rt_SIR_pcr.mat','s','t0','t1','q_mat','Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
data_real = load('results_Rt_SIR_real_pcr.mat','s','t0','t1','q_mat','Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
inputs = load('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_real','dI_inflow_smooth',...
    'pos_test_ratio_smooth','obs_ratio','obs_ratio_real','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1');

%% handle inputs
s = data.s; % s = setparam
t1 = dd(2020,11,30);
startFcast = t1+1;
endHist = startFcast-1;
startHist = data.t0;
startEstim = s.env_from;
startEstimFull = s.wave_2_from;
delay = 14;

%% handle mobility 
fcastPer = endFcast-startFcast+1;
mobilityFcast = forecast_mobility(mob,fcastPer,startHist,endHist,startWave);
% PCR data only	- as reported
Rt_data = tseries(data.t0+1:data.t1,data.q_mat(s.quant_idx_central,:));
Rt_data = resize(Rt_data,data.t0+1:t1);
% PCR data only	- as implied by cases in hospitals
Rt_data_real = tseries(data_real.t0+1:data_real.t1,data_real.q_mat(s.quant_idx_central,:));
Rt_data_real = resize(Rt_data_real,data_real.t0+1:t1);
Rt_series.default = Rt_data;
Rt_series.real = Rt_data_real;
[mobilityParams,mobilityParams_real] = estimate_mobility(mobilityFcast,Rt_series,startEstim,delay,startEstimFull);

%% saving results
save('mobility_forecast.mat','mobilityFcast','startHist','fcastPer','startWave');
save('mobility_estimation.mat','mobilityParams','mobilityParams_real','startEstim','delay','startEstimFull');

