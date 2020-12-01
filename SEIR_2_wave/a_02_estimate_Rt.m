%% initialization
% initialize;

%% setup
disp_from = dd(2020,8,1);
indiff = true; 
cut = 0;
dt = 1;

%% load data
load('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_pcr_adj','dI_inflow_pcr_adj_smooth',...
    'pos_test_ratio_smooth','obs_ratio_smooth','dI_inflow_smooth','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1','h_t00');
tt0 = t0+dt;
model_fnc = @estimate_Rt_SIR;
disp_to = t1-1;

%% calculations
s = setparam();
inputs_fnc = struct();
inputs_fnc.I0 = I0;
inputs_fnc.obs_ratio = [];
inputs_fnc.asymp_ratio = [];
inputs_fnc.z = double(resize(dI_inflow_smooth,t0:t1));
[Rt,q_mat,Yt,Rt_last,Rt_dist,Rt_rnd] = model_fnc(inputs_fnc,s,true,true,true);

inputs_fnc.z = double(resize(dI_inflow_pcr_smooth,t0:t1));
inputs_fnc.obs_ratio = double(resize(obs_ratio_smooth,t0:t1));
inputs_fnc.asymp_ratio = double(resize(asymp_ratio_smooth,t0:t1));
[Rt_pcr,q_mat_pcr,Yt_pcr,Rt_last_pcr,Rt_dist_pcr,Rt_rnd_pcr]  = model_fnc(inputs_fnc,s,true,true,true);

%% plotting stuff
figure('Name','Effective reproduction number, mean (PCR/PCR+AG)');
Rt_smooth_series_pcr = tseries(t0+1:t1,Rt_pcr);
Rt_smooth_series = tseries(t0+1:t1,Rt);
plot(resize(Rt_smooth_series,disp_from:t1),'linewidth',2);hold on;
plot(resize(Rt_smooth_series_pcr,disp_from:t1),'linewidth',2);hold on;
title('Mean Rt (smooth inputs)');
legend({'PCR + AG','PCR only'});
grid on;
% 
plot_fanchart(q_mat_pcr,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, PCR only)');
plot_fanchart(q_mat,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, PCR+AG)');

%% saving stuff
x.Rt_smooth = Rt_smooth_series;
x.Rt_smooth_pcr = Rt_smooth_series_pcr;
dbsave(x,'results.csv');

Rt = tseries(t0+1:t1,Rt);
save('results_Rt.mat','s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
save('results_Rt_pcr.mat','s','t0','t1','q_mat_pcr',...
    'Rt_pcr','Yt_pcr','Rt_last_pcr','Rt_dist_pcr','Rt_rnd_pcr');

%% optional: statistics
% model training
% if s.model_seir
%     tm0 = dd(2020,9,15);
%     tm1 = t1;
%     % tune theta
%     x = tseries(t0+1:t1,Yt.Iat./Yt.Iot);
%     x = tseries(t0+1:t1,Yt.Iaint./(Yt.Iaint+Yt.Isint));
%     x_tar = resize(asymp_ratio_smooth,t0:t1);
%     dx = 100*(x./x_tar-1);
%     figure; 
%     plot(resize(dx,tm0:tm1),'linewidth',1); grid on;
%     ylabel('%');
%     title('inflow I_{asymp}/I_{obs} (%)');
%     % tune obs_ratio
%     x = tseries(t0+1:t1,Yt.Ioint./(Yt.Ioint+Yt.Iuint));
%     x_tar = resize(obs_ratio_smooth,t0:t1);
%     dx = 100*(x./x_tar-1);
%     figure; 
%     plot(resize(dx,tm0:tm1),'linewidth',1); grid on;
%     ylabel('%');
%     title('inflow I_{obs}/I (%)');
% end