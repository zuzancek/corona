function [s] = setparam()

s = struct;
s.pop_size = 5443120;
s.sim_num = 10^4;
s.T = 100;
s.model_seir = true;
s.tests_min = 4700;
s.cases_min = 500;
s.ratio_threshold = 0.05;
s.env_from = dd(2020,10,8);
s.wave_2_from = dd(2020,9,1);

% serial interval (generation period)
s.SI.mean = 7.5;                s.SI.std = 0.62;
% time to test (observation period, from symptoms onset): "steady_state value"
s.T_test.mean = 3;             s.T_test.std = s.SI.std;      
% incubation period 
s.T_inc.mean = 5.2;             s.T_inc.std = s.SI.std;
% infectious period
s.T_inf.mean = 4.3;             s.T_inf.std = 0.62;
s.T_inf_asymp.mean = 4.3;       s.T_inf_asymp.std = 0.62;
s.T_inf_symp.mean = 4.3;        s.T_inf_symp.std = 0.62;
s.T_inf_obs.mean = 4.3;         s.T_inf_obs.std = 0.62;
s.T_inf_unobs.mean = 4.3;       s.T_inf_unobs.std = 0.62;
% sickness/symptoms period
s.T_sick_y.mean = 6.5;            s.T_sick_y.std = 0.62;
s.T_sick_o.mean = 9.5;            s.T_sick_o.std = 0.62;
% presymptomatic period 
s.T_pre.mean = s.T_inc.mean+s.T_inf.mean-s.SI.mean;           
s.T_pre.std = s.SI.std;
s.T_inf_obs0.mean = s.T_inf_obs.mean-s.T_pre.mean-s.T_test.mean;
s.T_inf_obs0.std = s.T_inf_obs.std;
s.SI_obs.mean = s.SI.mean-s.T_inf_obs.mean;
s.SI_obs.std = 0.62;
% latent period
s.T_lat.mean = s.T_inc.mean-s.T_pre.mean; s.T_lat.std  = s.SI.std;
s.share_reas = 1;
% days prior to hospital admission (from onset,adjusted for test period)
s.T_hosp.mean = 4.5;              s.T_hosp.std = 0.62;
s.T_hosp0.mean = s.T_hosp.mean-s.T_test.mean;    s.T_hosp0.std = 0.62;
% time to death (from ospital admission)
s.T_death.mean = 7;             s.T_death.std = 0.62;
% days at hospital (in case of recovery)
s.T_rec = 11.7;
% share of symptomatic patients in observed cases
s.symp_ratio_obs = 0.55;
% share of symptomatic patients needed to be hospitalized 
s.lambda = 0.0917; %0.0637; 
s.alpha_hr = 7.73/100; % recovery rate in hospital
s.eta_hr = s.alpha_hr.*s.T_rec;
s.alpha_ih = 6.37/100*(1/s.T_hosp.mean);
s.alpha_ir = 1/s.T_inf.mean-s.alpha_ih;
s.alpha_hd = 1.69/100; %1.23/100; % death rate in hospital
s.p_a_s = s.symp_ratio_obs*s.T_inf_asymp.mean/(s.T_inf_symp.mean);
% ICU rate
s.iota = 0.1375;
% discharge rate (non-ventilation case)
s.gamma_novent = 1/10; 
% case fatality rate (non-ventilation)
s.omega_novent = 0.1102; % 1-0.8125
% time from hospital (non-vent) admission to death
s.psi_novent = 1/6; 
% patients with ventilation needed (inflow)
s.xi = s.iota*0.5;
% time from vent-admission to death
s.psi_vent = 1/4;
s.omega_vent = 0.33;
% death probability/time
s.omega_y = 5.15/100;       s.omega_o = (37.16/100);  %alpha_hdy/alpha_hdo = 9.93 pct
s.T_death_y = 6.05;         s.T_death_o = 4.42;  
% recovery at hospital
s.T_rec_y = 9.2;            s.T_rec_o = 14.5;
% hospitalization probability/time to
s.eta_y = 2.32/100;         s.eta_o = 31.86/100;
s.eta_y = 2.94/100;         s.eta_o = 26.31/100;
s.T_hosp_y = 5;             s.T_hosp_o = 5;
s.T_hosp_y = 7;             s.T_hosp_o = 4;
% s.eta_y = 4.8/100;          s.eta_o = 15.6/100;
% s.eta_y = 5.42/100;         s.eta_o = 16.79/100;
s.alpha_weight = 0.25;
s.kappa_res_0 = 1/3;
s.kappa_res_delta_0 = -0.5;
s.kappa_res_alpha = 3;
s.beta_res = s.kappa_res_0/(1+s.kappa_res_delta_0^s.kappa_res_alpha);
s.kappa_res_0 = struct;
s.kappa_res_0.at = 0;
s.kappa_res_0.delta = 0;

g1.r0_scale = 10^2;
g1.mean = 0.98;
g1.alpha = 1;
s.g1 = g1;
g2.a0 = g1.mean+1/g1.r0_scale;
g2.a1 = 14;
g2.scale = 2.3;
g2.alpha = 1;
s.g2 = g2;
s.w_vec_default = 0.5+zeros(s.T,1);

s.old_share = 0.1385;
s.old_death_ratio = 0.8407;
s.obs_ratio_tar = 1/5;
s.obs_ratio = s.obs_ratio_tar;
s.self_isolation_effect = 1-0.05;
s.case_isolation_effect = 1/(1-1/3);
s.threshold = 0.05;
s.scale_fact = 4;

s.smooth_width = 7;
s.smooth_width_hosp = 5;
s.smooth_type = 5;
s.smooth_ends = 1;

s.quant = [0.05:0.05:0.95]; %#ok<*NBRAK>
s.quant_idx_central = 10;
s.quant_legend = cellstr(strcat(num2str((s.quant*100)'),'%'));

s.color_bkg = [0.176 0.176 .427];
s.color_grid = [0.69 0.91 0.973];
s.color_graph = [0.122 0.478 0.847];

pd = makedist('HalfNormal','mu',0,'sigma',1);
xx = 0:5; % last 6 days
weights = pdf(pd,xx);
s.pweight = weights/sum(weights);

% Rt distribution (kernel-smoothing & time-dependent)
s.min_pts = 100;
s.max_pts = 200;
s.min_dif = 1;
s.max_dif = 4;
s.shift_max = 2*s.SI.mean;

s.smoothing_method_data = @mov_median;
s.smoothing_method_params = @smooth_series;

end

