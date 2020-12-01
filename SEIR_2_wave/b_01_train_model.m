%% initialization
initialize;

%% load data
% inputs, ratios (asympt, obs), hospitals
data = load('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_pcr_adj',...
    'dI_inflow_smooth','dI_inflow_pcr_adj_smooth',...
    'pos_test_ratio_smooth','obs_ratio_smooth','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1','h_t00');
% reproduction number
data_Rt = load('results_Rt.mat','s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
data_Rt = load('results_Rt.mat','s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
% mobility
load('mobility_forecast.mat','mobilityFcast','startHist','fcastPer','startWave');
load('mobility_estimation.mat','mobilityParams','startEstim','delay','startEstimFull');

%% training dates
train_from = dd(2020,3,19);
train_to = dd(2020,7,1);
time_interval.dateFrom = train_from;
time_interval.dateTo = train_to;

%% handle inputs
% ******* Mobility
% training: included in effective Rt, set to 1
% forecast: medium scenario assumed
% mobility.values = mf.mobilityFcast.medium(startFcast:end)/100;
% mobility.x_grid = me.mobilityParams.f.x_grid;
% mobility.y_grid = me.mobilityParams.f.y_grid;
% mobility.scale = mf.mobilityFcast.medium(startFcast);
mobility.forecast = false;
% ******* Restrictions 
% training: included in effective Rt, set to 1
% forecast: unchanged after dateFrom (NPC assumed)
% restrictions = s.kappa_res_0; % delta, at;
restrictions.forecast = false;
% ******* Testing
obs_ratio = double(resize(data.obs_ratio_smooth,train_from:train_to));
asymp_ratio = double(resize(data.asymp_ratio_smooth,train_from:train_to));
T_test = 1+zeros(train_to-train_from+1,1);
% ******* Initial values
init.S = tseries(data.t0:data.t1-1,data_Rt.Yt.St);
init.Io = tseries(data.t0:data.t1-1,data_Rt.Yt.Iot);
init.I = tseries(data.t0:data.t1-1,data_Rt.Yt.It);
% ******* collect
inputs.init = init;
inputs.T_test = T_test;
inputs.obs_ratio = obs_ratio;
inputs.asymp_ratio = asymp_ratio;
inputs.mobility = mobility;
inputs.restrictions = restrictions;
inputs.Rt = tseries(data.t0:data.t1-1,data_Rt.Rt_rnd);

%% model training
[res_mean,res_quant] = train_SEIR(time_interval,inputs,data.s);

%% display results
tm0 = dd(2020,9,15);
tm1 = t1;
% tune theta
x = tseries(t0:t1-1,Yt.Ist./Yt.Iat);
as = as(2:end)/100;
x_tar = as./(1-as);
dx = 100*(x./x_tar-1);
figure; 
plot(resize(dx,tm0:tm1),'linewidth',1); grid on;
ylabel('%');
title('I_{symp}/I_{asymp} (%)');
% tune obs_ratio
x = tseries(t0:t1-1,Yt.Iot./Yt.It);
x_tar = s.obs_ratio_tar;
dx = 100*(x./x_tar-1);
figure; 
plot(resize(dx,tm0:tm1),'linewidth',1); grid on;
ylabel('%');
title('I_{obs}/I (%)');