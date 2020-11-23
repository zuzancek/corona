function [s] = setparam()

s = struct;
s.pop_size = 5443120;
s.sim_num = 10^4;
s.T = 100;
s.model_seir = true;
s.tests_min = 4700;
s.ratio_threshold = 0.1;
s.env_from = dd(2020,10,8);
s.wave_2_from = dd(2020,9,1);

% serial interval
extra = 0;
s.SI.mean = 6.5+extra;      s.SI.std = 0.62;
% time to test (observation period, from symptoms onset): "steady_state value"
s.T_test0 = 1;        
% incubation period 
s.T_inc.mean = 5;         s.T_inc.std = s.SI.std;
s.T_inf.mean = 3+extra;   s.T_inf.std = 0.62;
% presymptomatic period 
s.T_pre.mean = s.T_inc.mean+s.T_inf.mean-s.SI.mean;           
s.T_pre.std = s.SI.std;
% latent period
s.T_lat.mean = s.T_inc.mean-s.T_pre.mean; s.T_lat.std  = s.SI.std;
% infectious period - observed asymptomatic cases % 3.5
s.T_inf_asymp.mean = 3+extra;   s.T_inf_asymp.std = 0.62;
% infectious period - observed symptomatic cases
s.T_inf_symp.mean = 3.5+extra;    s.T_inf_symp.std = 0.62;
% infectious period - unobserved cases 3.75
s.share_reas = 0.8;
% days prior to hospital admission (from onset,adjusted for test period)
s.T_inf_hosp.mean = 6-s.T_test0;    s.T_inf_hosp.std = 0.62;
s.zeta = 1/s.T_inf_hosp.mean;
% share of symptomatic patients in observed cases
s.symp_ratio_obs = 0.43;
% share of symptomatic patients needed to be hospitalized 
s.lambda = 0.037; %0.0743; % <-- test here higher rate
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

s.obs_ratio_tar = 1/10;
s.obs_ratio = s.obs_ratio_tar;
s.self_isolation_effect = 1-0.12;
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

end

