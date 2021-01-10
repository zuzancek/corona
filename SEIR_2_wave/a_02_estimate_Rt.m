%% initialization
% initialize;
rng(1000);

%% setup
disp_from = dd(2020,9,1);
indiff = true; 
cut = 0;
dt = 1;

%% load data
load('inputs.mat','dates','cases_data','hosp_data','deaths_data','mob_data','s');

t0 = dates.t0;
t1 = dates.t1;
t2 = dates.t2;
tt0 = t0+dt;
s = setparam();
s.model_seir = false;
if s.model_seir
    model_fnc = @estimate_Rt_SEIR;
else
    model_fnc = @estimate_Rt_SIR;
end
disp_to = t2-1;

%% inputs definition
inputs_fnc = struct();
init = struct();
init.I0 = cases_data.I0;
init.H0 = hosp_data.H_smooth(t1-1);
init.D0 = hosp_data.D_smooth(t1);
params = struct();
params.obs_ratio = cases_data.obs_ratio;           
params.old_ratio = cases_data.old_ratio_smooth;
params.old_ratio = smooth_series(params.old_ratio);
params.death_ratio = deaths_data.old_ratio_smooth;
params.asymp_ratio = cases_data.asymp_ratio_smooth;

%% calculations
s = setparam();
% _mm = moving median; _smooth = quasi-gaussian smoother
% reported data, PCR only, testing is optimal (sstate observ.ratio)
inputs_fnc.z = double(resize(cases_data.cases_pcr_smooth,t0:t2));
inputs_fnc.I0 = cases_data.I0;
inputs_fnc.H0 = 0;
inputs_fnc.D0 = 0;
%
inputs_fnc.obs_ratio = [];
inputs_fnc.old_ratio = [];
inputs_fnc.death_ratio = []; 
inputs_fnc.asymp_ratio = [];
[Rt_pcr,q_mat_pcr,Yt_pcr,Rt_last_pcr,Rt_dist_pcr,Rt_rnd_pcr]  = model_fnc(inputs_fnc,s,true,true,true);
% reported data, PCR only, realisting testing (realistic observ.ratio)
inputs_fnc.obs_ratio = double(resize(cases_data.obs_ratio,t0:t2));
inputs_fnc.old_ratio = params.old_ratio;
inputs_fnc.death_ratio = params.death_ratio;
inputs_fnc.asymp_ratio = double(resize(cases_data.asymp_ratio,t0:t2));
[Rt_test,q_mat_test,Yt_test,Rt_last_test,Rt_dist_test,Rt_rnd_test]  = model_fnc(inputs_fnc,s,true,true,true);
% real data, PCR only, optimal testing
inputs_fnc.z = double(resize(cases_data.cases_pcr_implied,t0:t2));
inputs_fnc.obs_ratio = [];
inputs_fnc.old_ratio = params.old_ratio;
inputs_fnc.death_ratio = params.death_ratio;
inputs_fnc.asymp_ratio = [];
[Rt_real,q_mat_real,Yt_real,Rt_last_real,Rt_dist_real,Rt_rnd_real] = model_fnc(inputs_fnc,s,true,true,false);
% reported data, PCR+AG
inputs_fnc.z = double(resize(cases_data.cases_total_smooth,t0:t2));
inputs_fnc.obs_ratio = [];
inputs_fnc.old_ratio = [];
inputs_fnc.asymp_ratio = [];
inputs_fnc.death_ratio = [];
% I0 = Yt_pcr.Iot(t1-t0+1);
% inputs_fnc.I0 = I0;
[Rt_total,q_mat_total,Yt_total,Rt_last_total,Rt_dist_total,Rt_rnd_total]  = model_fnc(inputs_fnc,s,true,true,false);

%% plotting stuff
% 0./ cases
figure('Name','New cases (smooth data, means)')
plot(cases_data.cases_pcr_implied,'linewidth',1); hold on;
plot(cases_data.cases_pcr_smooth,'k--', 'linewidth',1);
plot(cases_data.cases_total_smooth,'Color',[0.5 0 0.5],'linewidth',1);
grid on;
title('New cases');
legend({'Implied by hospitals','PCR reported','AG+PCR reported'});

% 1./ reproduction number
figure('Name','Effective reproduction number, means');
nn = length(Rt_pcr);
Rt_smooth_series_pcr = tseries(t0+1:t2,Rt_pcr);
Rt_smooth_series_real = tseries(t0:t2,Rt_real);
Rt_smooth_series_test = tseries(t0+1:t2,Rt_test);
Rt_smooth_series_total = tseries(t0+1:t2,Rt_total);
plot(resize(Rt_smooth_series_pcr,disp_from:t2),'linewidth',2);hold on;
plot(resize(Rt_smooth_series_test,disp_from:t2),'linewidth',2);hold on;
plot(resize(Rt_smooth_series_real,disp_from:t2),'linewidth',2);hold on;
plot(resize(Rt_smooth_series_total,disp_from:t2),'linewidth',2);hold on;
title('Effective reproduction number (smooth inputs)');
legend({'reported data, PCR only, testing is optimal',...
    'reported data, PCR only, testing is realistic',...
    'implied data', 'reported data, AG+PCR'});
grid on;
% 
plot_fanchart(q_mat_real,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, implied data)',true);
plot_fanchart(q_mat_pcr,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, PCR only, reported data)',true);
plot_fanchart(q_mat_total,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, PCR+AG, reported data)',true);

% 2./ situation in hospitals
figure('Name','Hospitals & Deaths, means');
subplot(2,1,1)
Dt_pcr = tseries(t0+1:t2,Yt_pcr.Dt(1:nn));
Dt_real = tseries(t0+1:t2,Yt_real.Dt(1:nn));
Dt_total = tseries(t0+1:t2,Yt_total.Dt(1:nn));
plot(resize(Dt_pcr,disp_from:t2),'linewidth',2);hold on;
plot(resize(Dt_real,disp_from:t2),'linewidth',2);hold on;
plot(resize(Dt_total,disp_from:t2),'linewidth',2);hold on;
plot(resize(hosp_data.D_smooth,disp_from:t2),'k--','linewidth',1);hold on;
title('Deaths (smooth inputs)');
legend({'reported data, PCR only',...
    'implied data', 'reported data, PCR+AG'});
grid on;
subplot(2,1,2)
Ht_pcr = tseries(t0+1:t2,Yt_pcr.Ht(1:nn));
Ht_real = tseries(t0+1:t2,Yt_real.Ht(1:nn));
Ht_total = tseries(t0+1:t2,Yt_total.Ht(1:nn));
plot(resize(Ht_pcr,disp_from:t2),'linewidth',2);hold on;
plot(resize(Ht_real,disp_from:t2),'linewidth',2);hold on;
plot(resize(Ht_total,disp_from:t2),'linewidth',2);hold on;
plot(resize(hosp_data.H_smooth,disp_from:t2),'k--','linewidth',1);hold on;
title('Hospitals (smooth inputs)');
legend({'reported data, PCR only',...
    'implied data', 'reported data, PCR+AG'});
grid on;

%% saving stuff
x.Rt = Rt_smooth_series_pcr;
x.Rt_real = Rt_smooth_series_real;
x.Rt_test = Rt_smooth_series_test;
dbsave(x,'results.csv');

Rt = tseries(t0:t2-1,Rt_pcr); %#ok<*NASGU>
q_mat = q_mat_pcr; Yt = Yt_pcr; Rt_last = Rt_last_pcr; Rt_dist = Rt_dist_pcr; Rt_rnd = Rt_rnd_pcr;
save(strcat('results_Rt.mat'),'s','t0','t1','t2','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
Rt = tseries(t0:t2,Rt_real);
q_mat = q_mat_real; Rt = Rt_real; Yt = Yt_real; Rt_last = Rt_last_real; Rt_dist = Rt_dist_real; Rt_rnd = Rt_rnd_real;
save(strcat('results_Rt_real.mat'),'s','t0','t1','t2','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
Rt = tseries(t0:t2-1,Rt_total);
q_mat = q_mat_total; Yt = Yt_total; Rt_last = Rt_last_total; Rt_dist = Rt_dist_total; Rt_rnd = Rt_rnd_total;
save(strcat('results_Rt_total.mat'),'s','t0','t1','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
