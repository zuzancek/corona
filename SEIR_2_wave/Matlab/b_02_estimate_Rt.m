%% initialization
initialize;

%% setup
disp_from = dd(2020,9,1);
indiff = true;
cut = 0;
dt = 1;

%% load data
load('results/inputs.mat','dates','cases_data','hosp_data','deaths_data','mob_data','s');
load('results/results_impl.mat','cases_implied_data','hosp_implied_counterfact_data','implied_data','hosp_implied_data','s');

t0 = dates.t0;
t1 = dates.t1;
t2 = dates.t2;
tt0 = t0+dt;
disp_to = t2;
s = setparam();
s.estimate_Rt = false;

if s.estimate_Rt
    model_fnc = @estimate_Rt_SIR;
    
    
    %% inputs definition
    inputs_fnc = struct();
    init = struct();
    init.I0 = cases_data.I0;
    params = struct();
    % params.obs_ratio = cases_implied_data.obs_ratio;
    params.old_ratio = cases_data.old_ratio;
    
    % d0 = dd(2020,09,7); d1 = dd(2020,10,17);
    %% calculations
    s = setparam();
    
    % 1. reported data, PCR only, testing is optimal
    cases_data.cases_pcr_smooth = smooth_series(cases_data.cases_pcr_smooth);
    inputs_fnc.z = double(resize(cases_data.cases_pcr_smooth,t0:t2));
    inputs_fnc.I0 = cases_data.I0;
    inputs_fnc.old_ratio = cases_data.old_ratio;
    inputs_fnc.eta_o = [];
    inputs_fnc.eta_y = [];
    [Rt_pcr,q_mat_pcr,Yt_pcr,Rt_last_pcr,Rt_dist_pcr,Rt_rnd_pcr]  = model_fnc(inputs_fnc,s,true,true,false);
    
    % real (implied) data
    z0 = resize(cases_data.cases_pcr_smooth,t0:t2);
    z0(end-length(cases_implied_data.X_smooth_all)+1:end) = cases_implied_data.X_smooth_all;
    inputs_fnc.z = double(z0);
    inputs_fnc.old_ratio = cases_implied_data.rho_real;
    inputs_fnc.eta_o = hosp_implied_data.eta_o;
    inputs_fnc.eta_y = hosp_implied_data.eta_y;
    [Rt_real,q_mat_real,Yt_real,Rt_last_real,Rt_dist_real,Rt_rnd_real] = model_fnc(inputs_fnc,s,true,true,false);
    
    cases.reported = cases_data.cases_pcr;
    cases.implied = cases_implied_data.X_smooth_all;
    
    nn = length(Rt_pcr);
    Rt_smooth_series_pcr = tseries(t0:t0+nn-1,Rt_pcr);
    Rt_smooth_series_real = tseries(t0:t0+nn-1,Rt_real);
    Rt_implied.mean = Rt_real;
    Rt_reported.mean = Rt_pcr;
    
else
    % load results from external files produced by R-codes
    % load data series containing mean&CI's for Rt calculated from
    % official/implied cases
    info = load_fanchart_tseries('src_dir', '../R/Rt_estimation/results',...
        'src_filenames', {'output_R_reported.csv','output_R_implied.csv'},...
        'tar_dir','results','tar_filenames', {'Rt_reported.csv','Rt_implied.csv'});
    
    cases.reported = info{1}.X_ts;
    cases.implied = info{2}.X_ts;
    
    Rt_reported.mean = info{1}.mean_ts;
    Rt_implied.mean = info{2}.mean_ts;
    
end

%% plotting stuff
% 0./ cases
figure('Name','New cases (smooth data, means)')
plot(resize(cases.implied,disp_from:t2),'linewidth',1); hold on;
plot(resize(cases.reported,disp_from:t2),'k--', 'linewidth',1);
grid on;
title('New cases');
legend({'Implied by hospitals','PCR reported','AG+PCR reported'});

% 1./ reproduction number
figure('Name','Effective reproduction number, means');
plot(resize(Rt_reported.mean,disp_from+30:t2),'linewidth',2);hold on;
plot(resize(Rt_implied.mean,disp_from+30:t2),'linewidth',2);hold on;
title('Effective reproduction number (smooth inputs)');
legend({'reported data (PCR)','implied data'});
grid on;

plot_fancharts_cmp(q_mat_pcr(5:15,:),q_mat_real(5:15,:),s,disp_from+30,disp_to,...
    'offsetdate',t0,'title','Effective reproduction numbers (reported vs.implied, 75% CI)','legend',{'reported data (PCR)','implied data'});
%
plot_fanchart(q_mat_real,s,dt,disp_from+30,disp_to,t0,'Effective reproduction number (Rt, implied data)',true);
plot_fanchart(q_mat_pcr,s,dt,disp_from+30,disp_to,t0,'Effective reproduction number (Rt, PCR only, reported data)',true);

%% saving stuff
x.Rt = Rt_smooth_series_pcr;
x.Rt_real = Rt_smooth_series_real;
% x.Rt_test = Rt_smooth_series_test;
dbsave(x,'results/results_Rt_all.csv');

Rt = tseries(t0+1:t2,Rt_pcr); %#ok<*NASGU>
q_mat = q_mat_pcr; Yt = Yt_pcr; Rt_last = Rt_last_pcr; Rt_dist = Rt_dist_pcr; Rt_rnd = Rt_rnd_pcr;
save(strcat('results/results_Rt.mat'),'s','t0','t1','t2','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
Rt = tseries(t0+1:t2,Rt_real);
q_mat = q_mat_real; Rt = Rt_real; Yt = Yt_real; Rt_last = Rt_last_real; Rt_dist = Rt_dist_real; Rt_rnd = Rt_rnd_real;
save(strcat('results/results_Rt_real.mat'),'s','t0','t1','t2','q_mat',...
    'Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
