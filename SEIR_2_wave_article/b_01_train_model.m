%% initialization
initialize;

%% load data
% inputs, ratios (asympt, obs), hospitals
data = load('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_real','dI_inflow_smooth',...
    'pos_test_ratio_smooth','obs_ratio','obs_ratio_real','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1');
% reproduction number (PCR tests exclusively, testing is optimal (sstate))
data_Rt = load('results_Rt_SIR_pcr.mat','s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
% reproduction number (PCR tests exclusively, testing is relistic (sstate))
data_Rt_real = load('results_Rt_SIR_real_pcr.mat','s','t0','t1','q_mat',...
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
dt = 1;

s = data.s;

%% handle inputs
% ******* Mobility
% training: included in effective Rt, set to 1
% forecast: medium scenario assumed
% mobility.values = mf.mobilityFcast.medium(startFcast:end)/100;
% mobility.x_grid = me.mobilityParams.f.x_grid;
% mobility.y_grid = me.mobilityParams.f.y_grid;
% mobility.scale = mf.mobilityFcast.medium(startFcast);
mobility.forecast = false;
mobility_real = mobility;
% ******* Restrictions 
% training: included in effective Rt, set to 1
% forecast: unchanged after dateFrom (NPC assumed)
% restrictions = s.kappa_res_0; % delta, at;
restrictions.forecast = false;
restrictions_real = restrictions;
% ******* Testing
obs_ratio = double(resize(data.obs_ratio,train_from_date:train_to_date));
obs_ratio_real = double(resize(data.obs_ratio_real,train_from_date:train_to_date));
asymp_ratio = (1-data.s.symp_ratio_obs)+0*double(resize(data.asymp_ratio_smooth,train_from_date:train_to_date));
asymp_ratio_real = double(resize(data.asymp_ratio_smooth,train_from_date:train_to_date));
T_test = 1+zeros(train_to-train_from+1,1);
% ******* Initial values
init.St = data_Rt.Yt.St(train_from:train_to);
init.Iot = data_Rt.Yt.Iot(train_from:train_to);
init.Iut = data_Rt.Yt.Iut(train_from:train_to);
init.It = data_Rt.Yt.It(train_from:train_to);
init.Iat = data_Rt.Yt.Iat(train_from:train_to);
init.Ist = data_Rt.Yt.Ist(train_from:train_to);
init.Ht = data_Rt.Yt.Ht(train_from:train_to);
init.Dt = data_Rt.Yt.Dt(train_from:train_to);
init_real.St = data_Rt_real.Yt.St(train_from:train_to);
init_real.Iot = data_Rt_real.Yt.Iot(train_from:train_to);
init_real.Iut = data_Rt_real.Yt.Iut(train_from:train_to);
init_real.It = data_Rt_real.Yt.It(train_from:train_to);
init_real.Iat = data_Rt_real.Yt.Iat(train_from:train_to);
init_real.Ist = data_Rt_real.Yt.Ist(train_from:train_to);
init_real.Ht = data_Rt_real.Yt.Ht(train_from:train_to);
init_real.Dt = data_Rt_real.Yt.Dt(train_from:train_to);
% ******* collect
tinputs.init = init;
tinputs.Rt = data_Rt.Rt_rnd;%(:,train_from-2*data.s.SI.mean:train_to-2*data.s.SI.mean);
tinputs_real.init = init_real;
tinputs_real.Rt = data_Rt_real.Rt_rnd;%(:,train_from-2*data.s.SI.mean:train_to-2*data.s.SI.mean);
assumptions.T_test = T_test;
assumptions.t_hosp = 59;
assumptions.obs_ratio = obs_ratio;
assumptions.obs_ratio_effect = (data.s.obs_ratio./obs_ratio-1)+1;
assumptions.asymp_ratio = asymp_ratio;
assumptions.mobility = mobility;
assumptions.restrictions = restrictions;
assumptions_real.T_test = T_test;
assumptions_real.t_hosp = 59;
assumptions_real.obs_ratio = obs_ratio_real;
assumptions_real.obs_ratio_effect = 0*(data.s.obs_ratio./obs_ratio_real-1)+1;
assumptions_real.asymp_ratio = asymp_ratio_real;
assumptions_real.mobility = mobility_real;
assumptions_real.restrictions = restrictions_real;

%% model training
s.quant = 0.25:0.025:0.75;
[res_mean,res_quant] = train_SIR(time_interval,assumptions,tinputs,data.s);
[res_mean_real,res_quant_real] = train_SIR(time_interval,assumptions_real,tinputs_real,data.s);

%% display results
% infectious in distributions
ydata.data = init.It; ydata.col = 'r'; ydata.leg = 'Data-implied';
plot_fanchart(res_quant.It,s,dt,train_from,train_to,train_from,'Total active infections (SIR, PCR only, testing is optimal)',true,1,ydata);
ydata.data = init.Iot; ydata.col = 'r'; ydata.leg = 'Data-implied';
plot_fanchart(res_quant.Iot,s,dt,train_from,train_to,train_from,'Observed active infections (SIR, PCR only, testing is optimal)',true,1,ydata);
ydata_real.data = init_real.It; ydata_real.col = 'r'; ydata_real.leg = 'Data-implied';
plot_fanchart(res_quant_real.It,s,dt,train_from,train_to,train_from,'Total active infections (SIR, PCR only, testing is realistic)',true,1,ydata);
ydata_real.data = init_real.Iot; ydata_real.col = 'r'; ydata_real.leg = 'Data-implied';
plot_fanchart(res_quant_real.Iot,s,dt,train_from,train_to,train_from,'Observed active infections (SIR, PCR only, testing is realistic)',true,1,ydata);

% Data vs model-implied infectious
figure('Name','Data vs. model-implied observed infectious');
subplot(2,1,1);
plot(init.Iot,'linewidth',1);hold on;
plot(res_mean.Iot,'linewidth',1);hold on;
plot(res_quant.Iot(ceil(length(s.quant)/2),:),'linewidth',1);hold on;
grid on;
legend({'data implied','model-implied mean','model-implied median'});
title('Observed Active Infections (Optimal testing assumption)');

subplot(2,1,2);
plot(init_real.Iot,'linewidth',1);hold on;
plot(res_mean_real.Iot,'linewidth',1);hold on;
plot(res_quant_real.Iot(ceil(length(s.quant)/2),:),'linewidth',1);hold on;
grid on;
legend({'data implied','model-implied mean','model-implied median'});
title('Observed Active Infections (realistic testing assumption)');

%% save results
means = res_mean; init = tinputs; quants = res_quant; %#ok<*NASGU>
save('SIR.mat','means','quants','time_interval','assumptions','init');
means = res_mean_real; init = tinputs_real;quants = res_quant_real; assumptions = assumptions_real;
save('SIR_real.mat','means','time_interval','assumptions','init');
