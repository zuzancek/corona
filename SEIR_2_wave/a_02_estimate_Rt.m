%% initialization
initialize;

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
s = setparam();
% s.model_seir = false;
if s.model_seir
    model_fnc = @estimate_Rt_SEIR;
    del = 2;
else
    model_fnc = @estimate_Rt_SIR;
    del = 0;
end
disp_to = t1-del-1;

%% calculations
s = setparam();
inputs_fnc = struct();
inputs_fnc.I0 = I0;
inputs_fnc.sim_num = 1;
inputs_fnc.obs_ratio = double(resize(obs_ratio_smooth,t0:t1));
inputs_fnc.asymp_ratio = double(resize(asymp_ratio_smooth,t0:t1));
inputs_fnc.z = double(resize(dI_inflow_smooth,t0:t1));
[Rt_smooth,q_mat,Yt,x_mat,Rt_last] = model_fnc(inputs_fnc,s,true,true);

inputs_fnc.z = double(resize(dI_inflow_smooth,t0:t1));
inputs_fnc.obs_ratio = [];
inputs_fnc.asymp_ratio = [];
[Rt_adj_smooth,q_mat_adj,~,x_mat_adj] = model_fnc(inputs_fnc,s,true,true);

%% plotting stuff
figure;
Rt_adj_series = tseries(t0+1:t1-del,Rt_adj_smooth);
Rt_smooth_series = tseries(t0+1:t1-del,Rt_smooth);
plot(resize(Rt_smooth_series,disp_from:t1-del),'linewidth',1);hold on;
plot(resize(Rt_adj_series,disp_from:t1-del),'linewidth',1);hold on;
title('Rt (smooth inputs)');
legend({'observed inputs','adjusted'});
grid on;
% 
plot_fanchart(q_mat,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt)');
plot_fanchart(q_mat_adj,s,dt,disp_from,disp_to,t0,'Adjusted effective reproduction number (Rt)');
plot_fanchart(x_mat,s,dt,disp_from,disp_to,t0,'Active infections estimate (unobs.included)');
plot_fanchart(x_mat_adj,s,dt,disp_from,disp_to,t0,'Adjusted active infections estimate (unobs.included)');

%% saving stuff
x.Rt_smooth = Rt_smooth_series;
x.Rt_adj_smooth = Rt_adj_series;
dbsave(x,'results.csv');

Rt = tseries(t0+1:t1-del,Rt_smooth);
save('results_Rt.mat','q_mat','Rt','Yt','s','Rt_last','t0','t1');

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