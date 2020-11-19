%% initialization
initialize;

%% setup
disp_from = dd(2020,4,1);
indiff = true; 
cut = 0;
dt = 1;

%% load data
load('inputs.mat','dI_inflow','dI_inflow_smooth','dI_inflow_adj','dI_inflow_adj_smooth',...
    'pos_test_ratio','pos_test_ratio_smooth','I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1','h_t00');
as_data = dbload('data/asymp.csv');
% x = load('inputs.mat','dI_inflow','dI_inflow_smooth','dI_inflow_adj','dI_inflow_adj_smooth',...
%     'pos_test_ratio','pos_test_ratio_smooth','I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
%     'death_smooth','h_t0','h_t1','h_t00');
tt0 = t0+dt;
t1 = dd(2020,11,14);

model_fnc = @estimate_Rt_SEIR_aug;
del = 1;
disp_to = t1-del-1;

%% calculations
s = setparam();
% s.p_a_s = s.symp_ratio_obs*s.T_inf_asymp.mean/(s.T_inf_symp.mean*s.lambda+(1-s.lambda)*s.T_inf_hosp.mean);
% [Rt,~,~,Xt] = model_fnc(double(resize(dI_inflow_smooth,t0:t1)),I0,s);
[Rt_smooth,q_mat,Yt,x_mat,Rt_last] = model_fnc(double(resize(dI_inflow_smooth,t0:t1)),I0,s,s.quant,double(as_data.AS),s.pweight);


%% saving stuff
Rt_vec_raw = zeros(t1-t0+1-del,1);
Rt_vec_raw(dt+1:end) = Rt;
Rt_vec_raw = tseries(t0:t1-del,Rt_vec_raw);
x.Rt_raw = Rt_vec_raw;

Rt_vec_smooth = zeros(t1-t0+1-del,1);
Rt_vec_smooth(dt+1:end) = Rt_smooth;
Rt_vec_smooth = tseries(t0:t1-del,Rt_vec_smooth);
x.Rt_smooth = Rt_vec_smooth;

dbsave(x,'results.csv');
Rt = tseries(t0+1:t1-del,Rt_smooth);
save('results_Rt.mat','q_mat','Rt','Yt','s','Rt_last','t0','t1');

%% model training
tm0 = dd(2020,9,15);
tm1 = t1;
% tune theta
x = tseries(t0:t1-1,Yt.Ist./Yt.Iat);
x_tar = s.symp_ratio_obs/(1-s.symp_ratio_obs);
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