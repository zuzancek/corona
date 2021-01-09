%% initialization
initialize;
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
params0 = struct();
params0.obs_ratio = cases_data.obs_ratio;           
params0.old_ratio = cases_data.old_ratio_smooth;
params0.asymp_ratio = cases_data.asymp_ratio_smooth;
params = struct();
params.obs_ratio = resize(cases_data.obs_ratio,t1-1:t2);           
params.old_ratio = resize(cases_data.old_ratio_smooth,t1-1:t2);
params.asymp_ratio = resize(cases_data.asymp_ratio_smooth,t1-1:t2);
inputs_fnc0.init = init;
inputs_fnc.init = init;
inputs_fnc0.params = params0;
inputs_fnc.params = params;

%% calculations
% s = setparam();
% _mm = moving median; _smooth = quasi-gaussian smoother
% reported data, PCR only, testing is optimal (sstate observ.ratio)
inputs_fnc.z = double(resize(cases_data.cases_pcr_smooth,t0:t2));
inputs_fnc.I0 = cases_data.I0;
inputs_fnc.obs_ratio = []; 
inputs_fnc.asymp_ratio = [];
[Rt_pcr,q_mat_pcr,Yt_pcr,Rt_last_pcr,Rt_dist_pcr,Rt_rnd_pcr]  = model_fnc(inputs_fnc,s,true,true,true);
% reported data, PCR only, realisting testing (realistic observ.ratio)
inputs_fnc.obs_ratio = double(resize(cases_data.obs_ratio,t0:t2));
inputs_fnc.asymp_ratio = double(resize(cases_data.asymp_ratio,t0:t2));
[Rt_test,q_mat_test,Yt_test,Rt_last_test,Rt_dist_test,Rt_rnd_test]  = model_fnc(inputs_fnc,s,true,true,true);
% real data, PCR only, optimal testing
inputs_fnc.z = double(resize(cases_data.cases_pcr_implied,t0:t2));
inputs_fnc.obs_ratio = [];
inputs_fnc.asymp_ratio = [];
% I0 = Yt_pcr.Iot(t1-t0+1);
% inputs_fnc.I0 = I0;
[Rt_real,q_mat_real,Yt_real,Rt_last_real,Rt_dist_real,Rt_rnd_real]  = model_fnc(inputs_fnc,s,true,true,false);
% reported data, PCR+AG
inputs_fnc.z = double(resize(cases_data.cases_total_smooth,t0:t2));
inputs_fnc.obs_ratio = [];
inputs_fnc.asymp_ratio = [];
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
Dt_test = tseries(t0+1:t2,Yt_test.Dt(1:nn));
Dt_real = tseries(t0+1:t2,Yt_real.Dt(1:nn));
Dt_total = tseries(t0+1:t2,Yt_total.Dt(1:nn));
plot(resize(Dt_pcr,disp_from:t2),'linewidth',2);hold on;
plot(resize(Dt_test,disp_from:t2),'linewidth',2);hold on;
plot(resize(Dt_real,disp_from:t2),'linewidth',2);hold on;
plot(resize(Dt_total,disp_from:t2),'linewidth',2);hold on;
plot(resize(hosp_data.D_smooth,disp_from:t2),'k--','linewidth',1);hold on;
title('Deaths (smooth inputs)');
legend({'reported data, optimal testing, PCR only',...
    'reported data, realistic testing, PCR only', 'implied data', 'reported data, PCR+AG'});
grid on;
subplot(2,1,2)
Ht_pcr = tseries(t0+1:t2,Yt_pcr.Ht(1:nn));
Ht_test = tseries(t0+1:t2,Yt_test.Ht(1:nn));
Ht_real = tseries(t0+1:t2,Yt_real.Ht(1:nn));
Ht_total = tseries(t0+1:t2,Yt_total.Ht(1:nn));
plot(resize(Ht_pcr,disp_from:t2),'linewidth',2);hold on;
plot(resize(Ht_test,disp_from:t2),'linewidth',2);hold on;
plot(resize(Ht_real,disp_from:t2),'linewidth',2);hold on;
plot(resize(Ht_total,disp_from:t2),'linewidth',2);hold on;
plot(resize(hosp_data.H_smooth,disp_from:t2),'k--','linewidth',1);hold on;
title('Hospitals (smooth inputs)');
legend({'reported data, optimal testing, PCR only',...
    'reported data, realistic testing, PCR only', 'implied data', 'reported data, PCR+AG'});
grid on;

%% saving stuff
x.Rt = Rt_smooth_series_pcr;
x.Rt_real = Rt_smooth_series_real;
x.Rt_test = Rt_smooth_series_test;
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
