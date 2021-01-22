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
tshift = 6;

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

d0 = dd(2020,09,7); d1 = dd(2020,10,9);
%% calculations
s = setparam();
% _mm = moving median; _smooth = quasi-gaussian smoother
% reported data, PCR only, testing is optimal (sstate observ.ratio)
xx = resize(cases_data.cases_pcr_smooth,t0:t2);
xx(d0:d1) = NaN; xx = interp(xx,t0:t2);
cases_data.cases_pcr_smooth = smooth_series(xx);
inputs_fnc.z = double(resize(cases_data.cases_pcr_smooth,t0:t2));
inputs_fnc.I0 = cases_data.I0;
inputs_fnc.H0 = 0;
inputs_fnc.D0 = 0;
inputs_fnc.obs_ratio = [];
inputs_fnc.old_ratio = [];
inputs_fnc.death_ratio = []; 
inputs_fnc.asymp_ratio = [];
[Rt_pcr,q_mat_pcr,Yt_pcr,Rt_last_pcr,Rt_dist_pcr,Rt_rnd_pcr]  = model_fnc(inputs_fnc,s,true,true,false);
% real (implied) data
z0 = resize(cases_data.cases_pcr_smooth,t0:t2);
z0(t1+3:t2) = cases_data.cases_pcr_implied(t1+3:t2);
z0(d0:d1) = NaN; z0 = interp(z0,t0:t2);
cases_data.cases_pcr_implied = smooth_series(z0);
inputs_fnc.z = double(resize(cases_data.cases_pcr_implied,t0:t2));
inputs_fnc.obs_ratio = [];
inputs_fnc.old_ratio = params.old_ratio;
inputs_fnc.death_ratio = params.death_ratio;
inputs_fnc.asymp_ratio = [];
[Rt_real,q_mat_real,Yt_real,Rt_last_real,Rt_dist_real,Rt_rnd_real] = model_fnc(inputs_fnc,s,true,true,false);

%% plotting stuff
% 0./ cases
figure('Name','New cases (smooth data, means)')
plot(resize(cases_data.cases_pcr_implied,disp_from:t2),'linewidth',1); hold on;
plot(resize(cases_data.cases_pcr_smooth,disp_from:t2),'k--', 'linewidth',1);
grid on;
title('New cases');
legend({'Implied by hospitals','PCR reported','AG+PCR reported'});

% 1./ reproduction number
figure('Name','Effective reproduction number, means');
nn = length(Rt_pcr);
Rt_smooth_series_pcr = tseries(t0:t2,Rt_pcr);
Rt_smooth_series_real = tseries(t0:t2,Rt_real);
plot(resize(Rt_smooth_series_pcr,disp_from:t2),'linewidth',2);hold on;
plot(resize(Rt_smooth_series_real,disp_from:t2-tshift),'linewidth',2);hold on;
title('Effective reproduction number (smooth inputs)');
legend({'reported data (PCR)','implied data'});
grid on;
% 
plot_fanchart(q_mat_real,s,dt,disp_from,disp_to-tshift,t0,'Effective reproduction number (Rt, implied data)',true);
plot_fanchart(q_mat_pcr,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt, PCR only, reported data)',true);

% 2./ situation in hospitals
figure('Name','Hospitals & Deaths, means');
subplot(2,1,1)
Dt_pcr = tseries(t0:t2,Yt_pcr.Dt(1:nn));
Dt_real = tseries(t0:t2,Yt_real.Dt(1:nn));
% Dt_total = tseries(t0+1:t2,Yt_total.Dt(1:nn));
kappa_d = resize(Dt_real,disp_from:t2)./resize(hosp_data.D_smooth,disp_from:t2);
Dt_real = resize(Dt_real,disp_from:t2)./kappa_d;
Dt_pcr = resize(Dt_pcr,disp_from:t2)./kappa_d;
plot(resize(Dt_pcr,disp_from:t2),'linewidth',2);hold on;
plot(resize(Dt_real,disp_from:t2),'linewidth',2);hold on;
% plot(resize(Dt_total,disp_from:t2),'linewidth',2);hold on;
plot(resize(hosp_data.D_smooth,disp_from:t2),'k--','linewidth',1);hold on;
title('Deaths (smooth inputs)');
legend({'reported data, PCR only',...
    'implied data'});
grid on;
subplot(2,1,2)
Ht_pcr = tseries(t0:t2,Yt_pcr.Ht(1:nn));
Ht_real = tseries(t0:t2,Yt_real.Ht(1:nn));
kappa_h = resize(Ht_real,disp_from:t2)./resize(hosp_data.H_smooth,disp_from:t2);
Ht_real = resize(Ht_real,disp_from:t2)./kappa_h;
Ht_pcr = resize(Ht_pcr,disp_from:t2)./kappa_h;
% Ht_total = tseries(t0+1:t2,Yt_total.Ht(1:nn));
plot(resize(Ht_pcr,disp_from:t2),'linewidth',2);hold on;
plot(resize(Ht_real,disp_from:t2),'linewidth',2);hold on;
% plot(resize(Ht_total,disp_from:t2),'linewidth',2);hold on;
plot(resize(hosp_data.H_smooth,disp_from:t2),'k--','linewidth',1);hold on;
title('Hospitals (smooth inputs)');
legend({'reported data, PCR only',...
    'implied data'});
grid on;

%% saving stuff
x.Rt = Rt_smooth_series_pcr;
x.Rt_real = Rt_smooth_series_real;
% x.Rt_test = Rt_smooth_series_test;
dbsave(x,'results.csv');

Rt = tseries(t0:t2-1,Rt_pcr); %#ok<*NASGU>
q_mat = q_mat_pcr; Yt = Yt_pcr; Rt_last = Rt_last_pcr; Rt_dist = Rt_dist_pcr; Rt_rnd = Rt_rnd_pcr;
save(strcat('results_Rt.mat'),'s','t0','t1','t2','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
Rt = tseries(t0:t2-1,Rt_real);
q_mat = q_mat_real; Rt = Rt_real; Yt = Yt_real; Rt_last = Rt_last_real; Rt_dist = Rt_dist_real; Rt_rnd = Rt_rnd_real;
save(strcat('results_Rt_real.mat'),'s','t0','t1','t2','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
