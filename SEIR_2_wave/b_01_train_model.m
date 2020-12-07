%% initialization
% initialize;

%% load data
% inputs, ratios (asympt, obs), hospitals
data = load('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_real','dI_inflow_smooth',...
    'pos_test_ratio_smooth','obs_ratio','obs_ratio_real','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1');
% reproduction number (PCR tests exclusively/both PCR and AG)
data_Rt = load('results_Rt_SIR_pcr.mat','s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
% mobility
load('mobility_forecast.mat','mobilityFcast','startHist','fcastPer','startWave');
load('mobility_estimation.mat','mobilityParams','startEstim','delay','startEstimFull');

%% training dates
train_from_date = dd(2020,9,1);
train_to_date = dd(2020,11,27-1);
train_from = train_from_date-data_Rt.t0+1;
train_to = train_to_date-data_Rt.t0+1;
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
obs_ratio = double(resize(data.obs_ratio_smooth,train_from_date:train_to_date));
asymp_ratio = double(resize(data.asymp_ratio_smooth,train_from_date:train_to_date));
T_test = 1+zeros(train_to-train_from+1,1);
% ******* Initial values
init.S = data_Rt.Yt.St(train_from:train_to);
init.Io = data_Rt.Yt.Iot(train_from:train_to);
init.I = data_Rt.Yt.It(train_from:train_to);
% ******* collect
inputs.init = init;
inputs.T_test = T_test;
inputs.obs_ratio = obs_ratio;
inputs.obs_ratio_effect = (data.s.obs_ratio./obs_ratio-1)+1;
inputs.asymp_ratio = asymp_ratio;
inputs.mobility = mobility;
inputs.restrictions = restrictions;
inputs.Rt = data_Rt.Rt_rnd;%(:,train_from-2*data.s.SI.mean:train_to-2*data.s.SI.mean);

%% model training
s.quant = 0.25:0.025:0.75;
[res_mean,res_quant] = train_SIR(time_interval,inputs,data.s);

%% display results
% infectious in distributions
ydata.data = init.I; ydata.col = 'r'; ydata.leg = 'Data-implied';
plot_fanchart(res_quant.It,s,dt,train_from,train_to,train_from,'Total active infections (SIR, PCR only)',true,8,ydata);
ydata.data = init.Io; ydata.col = 'r'; ydata.leg = 'Data-implied';
plot_fanchart(res_quant.Iot,s,dt,train_from,train_to,train_from,'Observed active infections (SIR, PCR only)',true,8,ydata);

% Data vs model-implied infectious
figure('Name','Data vs. model-implied observed infectious');
plot(init.Io,'linewidth',1);hold on;
plot(res_mean.Iot,'linewidth',1);hold on;
plot(res_quant.Iot(ceil(length(s.quant)/2),:),'linewidth',1);hold on;
grid on;
legend({'data implied','model-implied mean','model-implied median'});

% TODO
% komparacie: (SIR, epi)
% - observed/unobserved;
% - asymptomatical/symptomatical
% tm0 = dd(2020,9,15);
% tm1 = t1;
% % tune theta
% x = tseries(t0:t1-1,Yt.Ist./Yt.Iat);
% as = as(2:end)/100;
% x_tar = as./(1-as);
% dx = 100*(x./x_tar-1);
% figure; 
% plot(resize(dx,tm0:tm1),'linewidth',1); grid on;
% ylabel('%');
% title('I_{symp}/I_{asymp} (%)');
% % tune obs_ratio
% x = tseries(t0:t1-1,Yt.Iot./Yt.It);
% x_tar = s.obs_ratio_tar;
% dx = 100*(x./x_tar-1);
% figure; 
% plot(resize(dx,tm0:tm1),'linewidth',1); grid on;
% ylabel('%');
% title('I_{obs}/I (%)');