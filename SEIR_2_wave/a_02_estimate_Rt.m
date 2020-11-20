%% initialization
initialize;

%% setup
disp_from = dd(2020,4,1);
indiff = true; 
cut = 0;
dt = 1;

%% load data
load('inputs.mat','dI_inflow','dI_inflow_smooth','dI_inflow_adj','dI_inflow_adj_smooth',...
    'pos_test_ratio','pos_test_ratio_smooth','obs_ratio_smooth','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1','h_t00');
tt0 = t0+dt;
t1 = dd(2020,11,14);

if s.model_seir
    model_fnc = @estimate_Rt_SEIR_aug;
    del = 1;
else
    model_fnc = @estimate_Rt_SIR;
    del = 0;
end
disp_to = t1-del-1;

%% calculations
s = setparam();
inputs_fnc = struct();
inputs_fnc.I0 = I0;
inputs_fnc.z = double(resize(dI_inflow_smooth,t0:t1));
inputs_fnc.obs_ratio = obs_ratio_smooth;
inputs_fnc.asymp_ratio = asymp_ratio_smooth;
[Rt,~,~,Xt] = model_fnc(inputs_fnc,s);
[Rt_smooth,q_mat,Yt,x_mat,Rt_last] = model_fnc(inputs_fnc,s,s.quant,s.pweight);

inputs_fnc.z = double(resize(dI_inflow_smooth,t0:t1));
inputs_fnc.obs_ratio = obs_ratio_smooth;
[Rt_adj,~,~,Xt_adj] = model_fnc(inputs_fnc,s);
[Rt_adj_smooth,q_mat_adj,~,x_mat_adj] = model_fnc(inputs_fnc,s,s.quant);

pos_test_ratio_smooth = 0*pos_test_ratio+pos_test_ratio_smooth;

%% plotting stuff
y = 0*resize(dI_inflow,t0:t1-del);
figure;
subplot(2,1,1);
plot(y+Rt,'linewidth',1);hold on;
plot(y+Rt_smooth,'linewidth',1);hold on;
title('Rt');
legend({'raw','smooth'});
grid on;
subplot(2,1,2);
plot(y+Rt_adj,'linewidth',1);hold on;
plot(y+Rt_adj_smooth,'linewidth',1);hold on;
title('Rt (with assumption)');
legend({'raw','smooth'});
grid on;
figure;
plot(y+Rt_smooth,'linewidth',1);hold on;
plot(y+Rt_adj_smooth,'linewidth',1);hold on;
title('Rt smooth');
legend({'observed','adjusted'});
grid on;
% 
plot_fanchart(q_mat,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt)');
plot_fanchart(q_mat_adj,s,dt,disp_from,disp_to,t0,'Adjusted effective reproduction number (Rt)');
plot_fanchart(x_mat,s,dt,disp_from,disp_to,t0,'Active infections estimate (unobs.included)');
plot_fanchart(x_mat_adj,s,dt,disp_from,disp_to,t0,'Adjusted active infections estimate (unobs.included)');

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

%% optional: statistics
% model training
if s.model_seir
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
end