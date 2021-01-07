%% initialization
initialize;

%% setup
disp_from = dd(2020,9,1);
indiff = true; 
cut = 0;
dt = 1;

%% load data
load('inputs.mat','dates','cases_data','hosp_data','deaths_data','mob_data','s');

t0 = dates.t0;
t1 = dates.t1;
tt0 = t0+dt;
s.model_seir = true;
if s.model_seir
    model_fnc = @estimate_Rt_SEIR;
else
    model_fnc = @estimate_Rt_SIR;
end
disp_to = t1-1;

%% inputs definition
inputs_fnc = struct();
init = struct();
init.I0 = cases_data.I0;
init.H0 = hosp_data.H_smooth(t0-1);
init.D0 = hosp_data.D_smooth(t0-1);
params = struct();
params.obs_ratio = cases_data.obs_ratio;
params.old_ratio = cases_data.old_ratio_smooth;
params.asymp_ratio = cases_data.asymp_ratio_smooth;
inputs_fnc.init = init;
inputs_fnc.params = params;

%% calculations
% s = setparam();
% _mm = moving median; _smooth = quasi-gaussian smoother
% 1.) reported data, PCR tests only
inputs_fnc.z = double(resize(cases_data.cases_pcr_mm,t0:t1)); % moving median
[Rt,q_mat,Yt,Rt_last,Rt_dist,Rt_rnd] = model_fnc(inputs_fnc,s,true,true,false); %#ok<*ASGLU>

nn = length(Rt);
% reported data, PCR only, testing is optimal (sstate observ.ratio)
inputs_fnc.z = double(resize(dI_inflow_pcr_smooth,t0:t1));
inputs_fnc.obs_ratio = []; %double(resize(obs_ratio,t0:t1));
inputs_fnc.asymp_ratio = [];
[Rt_pcr,q_mat_pcr,Yt_pcr,Rt_last_pcr,Rt_dist_pcr,Rt_rnd_pcr]  = model_fnc(inputs_fnc,s,true,true,true);
% reported data, PCR only, realisting testing (realistic observ.ratio)
inputs_fnc.z = double(resize(dI_inflow_pcr_smooth,t0:t1));
inputs_fnc.obs_ratio = double(resize(obs_ratio_real,t0:t1));
inputs_fnc.asymp_ratio = double(resize(asymp_ratio_smooth,t0:t1));
[Rt_test,q_mat_test,Yt_test,Rt_last_test,Rt_dist_test,Rt_rnd_test]  = model_fnc(inputs_fnc,s,true,true,true);
% real data, PCR only, optimal testing (it is correct!!!)
inputs_fnc.z = double(resize(dI_inflow_real,t0:t1));
inputs_fnc.obs_ratio = [];
inputs_fnc.asymp_ratio = double(resize(asymp_ratio_smooth,t0:t1));
[Rt_real,q_mat_real,Yt_real,Rt_last_real,Rt_dist_real,Rt_rnd_real]  = model_fnc(inputs_fnc,s,true,true,false);

%% plotting stuff
figure('Name','Effective reproduction number, means');
Rt_smooth_series_pcr = tseries(t0+1:t1,Rt_pcr(1:nn));
Rt_smooth_series_real = tseries(t0+1:t1,Rt_real(1:nn));
Rt_smooth_series_test = tseries(t0+1:t1,Rt_test(1:nn));
Rt_smooth_series = tseries(t0+1:t1,Rt(1:nn));
plot(resize(Rt_smooth_series_pcr,disp_from:t1),'linewidth',2);hold on;
plot(resize(Rt_smooth_series_test,disp_from:t1),'linewidth',2);hold on;
plot(resize(Rt_smooth_series,disp_from:t1),'linewidth',2);hold on;
title('Hospitalisations (smooth inputs)');
legend({'reported data, testing is optimal, PCR only',...
    'reported data, testing is realistic, PCR only',...
    'reported data, PCR+AG'});
grid on;
% 
figure('Name','Hospitals & Deaths, means');
subplot(2,1,1)
Dt_pcr = tseries(t0+1:t1,Yt_pcr.Dt(1:nn));
Dt_test = tseries(t0+1:t1,Yt_test.Dt(1:nn));
Ht_pcr = tseries(t0+1:t1,Yt_pcr.Ht(1:nn));
Ht_test = tseries(t0+1:t1,Yt_test.Ht(1:nn));
plot(resize(Dt_pcr,disp_from:t1),'linewidth',2);hold on;
plot(resize(Dt_test,disp_from:t1),'linewidth',2);hold on;
title('Deaths (smooth inputs)');
legend({'implied by reported data, testing is optimal, PCR only',...
    'implied by reported data, testing is realistic, PCR only'});
grid on;
subplot(2,1,2)
plot(resize(Ht_pcr,disp_from:t1),'linewidth',2);hold on;
plot(resize(Ht_test,disp_from:t1),'linewidth',2);hold on;
title('Hospitals (smooth inputs)');
legend({'implied by reported data, testing is optimal, PCR only',...
    'implied by reported data, testing is realistic, PCR only'});
grid on;
% 
plot_fanchart(q_mat_test,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, PCR only, realistic data/testing)',true);
plot_fanchart(q_mat_pcr,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, PCR only, as reported)',true);
plot_fanchart(q_mat,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, PCR+AG as reported)',true);

%% saving stuff
x.Rt_smooth = Rt_smooth_series;
x.Rt_smooth_pcr = Rt_smooth_series_pcr;
dbsave(x,'results.csv');

Rt = tseries(t0+1:t1,Rt); %#ok<*NASGU>
ext_opt = {'_SIR','_SEIR'};
save(strcat('results_Rt',ext_opt{1+s.model_seir},'.mat'),'s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
q_mat = q_mat_pcr; Rt = Rt_pcr; Yt = Yt_pcr; Rt_last = Rt_last_pcr; Rt_dist = Rt_dist_pcr; Rt_rnd = Rt_rnd_pcr;
save(strcat('results_Rt',ext_opt{1+s.model_seir},'_pcr.mat'),'s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
q_mat = q_mat_test; Rt = Rt_test; Yt = Yt_test; Rt_last = Rt_last_test; Rt_dist = Rt_dist_test; Rt_rnd = Rt_rnd_test;
save(strcat('results_Rt',ext_opt{1+s.model_seir},'_real_pcr.mat'),'s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
