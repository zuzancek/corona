%% initialization & cleanup
initialize;

%% setup
endFcast = dd(2020,12,31); % end of the year
startWave = dd(2020,09,01); % start of the second wave

%% loading stuff
x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility.csv','dateFormat','yyyy-mm-dd','freq','daily');
data = load('results_Rt.mat','q_mat_pcr','Rt_pcr','Yt','s','Rt_last','t0','t1');

%% handle inputs
s = data.s; % s = setparam
startFcast = data.t1+1;
endHist = startFcast-1;
startHist = data.t0;
startEstim = s.env_from;
startEstimFull = s.wave_2_from;
delay = 14;

%% handle mobility 
fcastPer = endFcast-startFcast+1;
mobilityFcast = forecast_mobility(mob,fcastPer,startHist,startWave);
% PCR data only	
Rt_data = tseries(data.t0+1:data.t1,data.q_mat_pcr(s.quant_idx_central,:));
mobilityParams = estimate_mobility(mobilityFcast,Rt_data,startEstim,delay,startEstimFull);

%% saving results
save('mobility_forecast.mat','mobilityFcast','startHist','fcastPer','startWave');
save('mobility_estimation.mat','mobilityParams','startEstim','delay','startEstimFull');

