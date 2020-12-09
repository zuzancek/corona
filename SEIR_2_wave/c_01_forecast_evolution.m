%% initialization & cleanup
initialize;

%% setup
endFcast = dd(2020,12,31); % end of the year
startWave = dd(2020,09,01); % start of the second wave

%% loading stuff
data = load('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_real','dI_inflow_smooth',...
    'pos_test_ratio_smooth','obs_ratio','obs_ratio_real','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1');
data_Rt = load('results_Rt_SIR_pcr.mat','s','t0','t1','q_mat','Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
data_real_Rt = load('results_Rt_SIR_real_pcr.mat','s','t0','t1','q_mat','Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
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
mobility.values = mf.mobilityFcast.medium(startFcast-mobility.delay-1:end-mobility.delay)/100;
mobility.scale = mf.mobilityFcast.medium(startFcast-mobility.delay)/100;
mobility.x_grid = me.mobilityParams.x;
mobility.y_grid_pos = me.mobilityParams.y_pos;
mobility.y_grid_neg = me.mobilityParams.y_neg;
% ****************** Restrictions
% unchanged for now
restrictions.forecast = false;
% ****************** Testing Quality
% current values/NPC
obs_ratio = data.obs_ratio(data.t1)+zeros(perFcast,1);
asymp_ratio = (1-data.s.symp_ratio_obs)+zeros(perFcast,1);
T_test = ones(perFcast,1);
% ******* Initial values
init.Iot = data_Rt.Yt.Iot(end);
init.Ht = data_Rt.Yt.Ht(end);
init.Dt = data_Rt.Yt.Dt(end);

%% forecast with SIRh
% put inputs together
inputs.init = init;
inputs.Rt = data_Rt.Rt_rnd;
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