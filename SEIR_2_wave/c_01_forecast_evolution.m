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
mf = load('mobility_forecast.mat','mobilityFcast','startHist','fcastPer','startWave');
me = load('mobility_estimation.mat','mobilityParams','mobilityParams_real','startEstim','delay','startEstimFull');

%% handle inputs
% ****************** Dates
s = data.s; % s = setparam
startFcast = data.t1+1;
endHist = startFcast-1;
startHist = data.t0;
startEstim = dd(2020,10,10);
startEstimFull = dd(2020,9,15);
endFcast = dd(2020,12,31);
perFcast = endFcast-startFcast+1;
% ****************** Mobility
% medium scenario assumed
mobility.forecast = true;
mobility.delay = me.mobilityParams.d;
mobility.values = mf.mobilityFcast.medium(startFcastmobility.delay:end-mobility.delay)/100;
mobility.scale = mf.mobilityFcast.medium(startFcast-mobility.delay);
mobility.x_grid = me.mobilityParams.x;
mobility.y_grid = me.mobilityParams.y_pos;
% ****************** Restrictions
% unchanged for now
restrictions.forecast = false;
% ****************** Testing Quality
% current values/NPC
obs_ratio = data.obs_ratio(data.t1)+zeros(perFcast,1);
asymp_ratio = (1-data.s.symp_ratio_obs)+zeros(perFcast,1);
T_test = ones(perFcast,1);
% ******* Initial values
init.Iot = data.Yt.Iot(data.t1);
init.Ht = data.Yt.Ht(data.t1);
init.Dt = data.Yt.Dt(data.t1);

%% forecast with SIRh
% put inputs together
inputs.init = init;
inputs.Rt = data.Rt_rnd;
assumptions.T_test = T_test;
assumptions.obs_ratio = obs_ratio;
assumptions.asymp_ratio = asymp_ratio;
assumptions.mobility = mobility;
assumptions.restrictions = restrictions;
assumptions_real.T_test = T_test;
time_interval.dateFrom = startFcast;
time_interval.dateTo = endFcast;

% run forecast
s.quant = 0.25:0.025:0.75;
[res_mean,res_quant] = train_SIR(time_interval,assumptions,inputs,data.s);

%% plot results
plot_fanchart(res_quant.dIt,s,1,startFcast,endFcast,startFcast);