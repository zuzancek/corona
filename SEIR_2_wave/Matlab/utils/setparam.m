function [s] = setparam(varargin)

s = struct;
s.pop_size = 5454147;
s.dep_ratio_65 = 16.6/100;
s.sim_num = 10^4;
s.T = 100;
s.estimate_Rt = false;
s.tests_min = 4700;
s.cases_min = 750;
s.ratio_threshold = 0.05;
s.env_from = dd(2020,10,8);
s.wave_2_from = dd(2020,9,1);
s.firstData_offset = 31;
s.deaths_with_covid_share = .18;
s.S_H_rate = 0.16;
s.S_H_rate_0 = 0.25;
s.scale_s_h = s.S_H_rate_0/s.S_H_rate;

% **** epidemiology
% serial interval (generation period)
s.SI.mean = 6.5;                s.SI.std = 0.62;
% incubation period
s.T_inc.mean = 5.1;             s.T_inc.std = s.SI.std;
s.obj_inc = makedist('Gamma','a',s.T_inc.mean*s.T_inc.std*s.T_inc.std,'b',1/(s.T_inc.std*s.T_inc.std));
% infectious period
s.T_inf.mean = 3.4;             s.T_inf.std = 0.62;
s.obj_inf = makedist('Gamma','a',s.T_inf.mean*s.T_inf.std*s.T_inf.std,'b',1/(s.T_inf.std*s.T_inf.std));
s.k_inf = 20;
s.pdf_inf = cut_tail(pdf(s.obj_inf,0:s.k_inf),5);
% presymptomatic period
s.k_pre = 15;
s.T_pre.mean = s.T_inc.mean+s.T_inf.mean-s.SI.mean;
s.T_pre.std = s.SI.std;
s.obj_pre = makedist('Gamma','a',s.T_pre.mean*s.T_pre.std*s.T_pre.std,'b',1/(s.T_pre.std*s.T_pre.std));
s.obj_pre_inf = makedist('Gamma','a',(s.T_pre.mean+s.T_inf.mean)*s.T_pre.std*s.T_pre.std,'b',1/(s.T_pre.std*s.T_pre.std));
% latent period
s.k_lat = 20;
s.T_lat.mean = s.T_inc.mean-s.T_pre.mean; s.T_lat.std  = s.SI.std;
s.obj_lat = makedist('Gamma','a',s.T_lat.mean*s.T_lat.std*s.T_lat.std,'b',1/(s.T_lat.std*s.T_lat.std));

% **** immunity, vaccination
s.psi_im_i = .25;  s.psi_im_h = .05;    s.psi_vac = .1;
s.phi_im_i = .25;  s.phi_im_h = .05;    s.phi_vac = .1;
s.T_im_i = 1.5*30; s.obj_im_i = makedist('Exponential','mu',s.T_im_i);
s.T_im_h = 4*30;   s.obj_im_h = makedist('Exponential','mu',s.T_im_h);
% **** mobility loss, quaranteen
s.alpha_s_o = 0.75; s.alpha_s_y = 1;
s.alpha_i_o = 0.25; s.alpha_i_y = 0.5;
s.mu = .85;

% **** testing
% time to test (observation period, from symptoms onset): "steady_state value"
s.T_test.mean = 2;              s.T_test.std = s.SI.std;
s.T_test_0 = 1;
s.k_test = 5;
s.k_pre_test = s.k_pre+s.k_test;
s.T_pre_test = s.T_test; s.T_pre_test.mean = s.T_pre_test.mean+s.T_pre.mean;
s.T_inf_obs.mean = 2;           s.T_inf_obs.std = s.SI.std; 
s.obj_pre_test = makedist('Gamma','a',s.T_pre_test.mean*s.T_pre_test.std*s.T_pre_test.std,'b',1/(s.T_pre_test.std*s.T_pre_test.std));
s.pdf_pre_test = cut_tail(pdf(s.obj_pre_test,0:s.k_inf),5);

% **** clinical characteristics
% sickness/symptoms period
s.T_sick_y = 8;                 s.T_sick_o = 8;       s.T_sick = 8;
s.T_sick_std = s.SI.std;
s.k_sick = 25;                  s.T_sick_pdf_type = 'Gamma';
try
    db = load('results/optimal_fit.mat','stat_total','stat_severe','stat_mild');
    set_prob_data();
catch err %#ok<NASGU>
    a_00_run_statistics();
    db = load('results/optimal_fit.mat','s');
    set_prob_data();
end

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

% big ratios (long-term averages)
s.symp_ratio_obs = 0.55;
s.old_share = 0.1385;
s.old_death_ratio = 0.8407;
s.obs_ratio_tar = 1/5;
s.obs_ratio = s.obs_ratio_tar;
s.self_isolation_effect = 1-0.05;
s.case_isolation_effect = 1/(1-1/3);
s.threshold = 0.05;
s.scale_fact = 4;

% smoothing
s.smooth_width = 7;
s.smooth_type = 5;
s.smooth_ends = 1;

% fancharts - quantiles
s.quant = [0.05:0.05:0.95]; %#ok<*NBRAK>
s.quant_idx_central = 10;
s.quant_legend = cellstr(strcat(num2str((s.quant*100)'),'%'));

s.quant_ext = [0.025 s.quant 0.975];
s.quant_ext_idx_central = 11;
s.quant_ext_legend = cellstr(strcat(num2str((s.quant*100)'),'%'));

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

s.smoothing_method_data = @mov_median_adj;
s.smoothing_method_params = @smooth_series;

    function[]=set_prob_data()
        % dbs definitions
        cutoff = 5;
        db_t = db.stat_total; db_s = db.stat_severe; db_m = db.stat_mild;
        % death (serious cases only, use db_s database)
        s.k_death = 40;
        s.omega_y = db_t.opt_fit_d_y.alpha;%/(1-s.deaths_with_covid_share);
        s.omega_o = db_t.opt_fit_d_o.alpha;%/(1-s.deaths_with_covid_share);
        s.pdf_hd_y = max(0,cut_tail(db_t.opt_fit_d_y.pdf(1:s.k_death+1),cutoff));
        s.pdf_hd_o = max(0,cut_tail(db_t.opt_fit_d_o.pdf(1:s.k_death+1),cutoff));
        s.epdf_hd_y = max(0,cut_tail(db_t.opt_fit_d_y.epdf(1:s.k_death+1),cutoff));
        s.epdf_hd_o = max(0,cut_tail(db_t.opt_fit_d_o.epdf(1:s.k_death+1),cutoff));
        s.T_death_y_mean = db_t.opt_fit_d_y.mean;
        s.T_death_o_mean = db_t.opt_fit_d_o.mean;
        s.time_d = reshape(db_t.opt_fit_d_y.time_grid(1:s.k_death+1),1,[]);
        s.omega_y_s = db_s.opt_fit_d_y.alpha/(1-s.deaths_with_covid_share);
        s.omega_o_s = db_s.opt_fit_d_o.alpha/(1-s.deaths_with_covid_share);
        s.pdf_sd_y = max(0,cut_tail(db_s.opt_fit_d_y.pdf(1:s.k_death+1),cutoff));
        s.pdf_sd_o = max(0,cut_tail(db_s.opt_fit_d_o.pdf(1:s.k_death+1),cutoff));
        s.epdf_sd_y = max(0,cut_tail(db_s.opt_fit_d_y.epdf(1:s.k_death+1),cutoff));
        s.epdf_sd_o = max(0,cut_tail(db_s.opt_fit_d_o.epdf(1:s.k_death+1),cutoff));
        s.obj_sd_y = db_s.opt_fit_d_y.eobj;
        s.obj_sd_o = db_s.opt_fit_d_o.eobj;
        s.obj_hd_y = db_t.opt_fit_d_y.eobj;
        s.obj_hd_o = db_t.opt_fit_d_o.eobj;
        s.T_death_y_mean_s = db_s.opt_fit_d_y.mean;
        s.T_death_o_mean_s = db_s.opt_fit_d_o.mean;
        s.time_d_s = reshape(db_s.opt_fit_d_y.time_grid(1:s.k_death+1),1,[]);
        % admission to ICU, ventilation, ECMO
        s.k_ser = 15;
        s.theta_y = 0.85*db_s.opt_fit_h_y.alpha;
        s.theta_o = 0.85*db_s.opt_fit_h_o.alpha;
        s.pdf_is_y = cut_tail(db_s.opt_fit_h_y.pdf(1:s.k_ser+1),cutoff);
        s.pdf_is_o = cut_tail(db_s.opt_fit_h_o.pdf(1:s.k_ser+1),cutoff);
        s.epdf_is_y = cut_tail(db_s.opt_fit_h_y.epdf(1:s.k_ser+1),cutoff);
        s.epdf_is_o = cut_tail(db_s.opt_fit_h_o.epdf(1:s.k_ser+1),cutoff);
        s.obj_is_o = db_s.opt_fit_h_y.eobj;
        s.obj_is_y = db_s.opt_fit_h_o.eobj;
        s.T_ser_y_mean = db_s.opt_fit_h_y.mean;
        s.T_ser_o_mean = db_s.opt_fit_h_o.mean;
        s.time_s = reshape(db_s.opt_fit_h_y.time_grid(1:s.k_ser+1),1,[]);
        % hospital admission
        s.k_hosp = 25;
        s.eta_y = 1.65*db_t.opt_fit_h_y.alpha; % 1.25
        s.eta_o = db_t.opt_fit_h_o.alpha; % 0.85
        s.pdf_ih_y = cut_tail(db_t.opt_fit_h_y.pdf(1:s.k_hosp+1),cutoff);
        s.pdf_ih_o = cut_tail(db_t.opt_fit_h_o.pdf(1:s.k_hosp+1),cutoff);
        s.epdf_ih_y = cut_tail(db_t.opt_fit_h_y.epdf(1:s.k_hosp+1),cutoff);
        s.epdf_ih_o = cut_tail(db_t.opt_fit_h_o.epdf(1:s.k_hosp+1),cutoff);
        s.obj_ih_y = db_t.opt_fit_h_y.eobj;
        s.obj_ih_o = db_t.opt_fit_h_o.eobj;
        s.T_hosp_y_mean = db_t.opt_fit_h_y.mean;
        s.T_hosp_o_mean = db_t.opt_fit_h_o.mean;
        s.time_h = reshape(db_t.opt_fit_h_y.time_grid(1:s.k_hosp+1),1,[]);
        % recovery (hospital (incl.serious),home)
        s.k_rec = 40;
        s.pdf_hr_y = cut_tail(db_t.opt_fit_r_y.pdf(1:s.k_rec+1),cutoff);
        s.pdf_hr_o = cut_tail(db_t.opt_fit_r_o.pdf(1:s.k_rec+1),cutoff);
        s.epdf_hr_y = cut_tail(db_t.opt_fit_r_y.epdf(1:s.k_rec+1),cutoff);
        s.epdf_hr_o = cut_tail(db_t.opt_fit_r_o.epdf(1:s.k_rec+1),cutoff);
        s.obj_hr_y = db_t.opt_fit_r_y.eobj;
        s.obj_hr_o = db_t.opt_fit_r_o.eobj;
        s.T_rec_h_y_mean = db_t.opt_fit_r_y.mean;
        s.T_rec_h_o_mean = db_t.opt_fit_r_o.mean;
        s.pdf_sr_y = cut_tail(db_s.opt_fit_r_y.pdf(1:s.k_rec+1),cutoff);
        s.pdf_sr_o = cut_tail(db_s.opt_fit_r_o.pdf(1:s.k_rec+1),cutoff);
        s.epdf_sr_y = cut_tail(db_s.opt_fit_r_y.epdf(1:s.k_rec+1),cutoff);
        s.epdf_sr_o = cut_tail(db_s.opt_fit_r_o.epdf(1:s.k_rec+1),cutoff);
        s.obj_sr_y = db_s.opt_fit_r_y.eobj;
        s.obj_sr_o = db_s.opt_fit_r_o.eobj;
        s.T_rec_s_y_mean = db_s.opt_fit_r_y.mean;
        s.T_rec_s_o_mean = db_s.opt_fit_r_o.mean;
        s.pdf_mr_y = cut_tail(db_m.opt_fit_r_y.pdf(1:s.k_rec+1),cutoff);
        s.pdf_mr_o = cut_tail(db_m.opt_fit_r_o.pdf(1:s.k_rec+1),cutoff);
        s.epdf_mr_y = cut_tail(db_m.opt_fit_r_y.epdf(1:s.k_rec+1),cutoff);
        s.epdf_mr_o = cut_tail(db_m.opt_fit_r_o.epdf(1:s.k_rec+1),cutoff);
        s.T_rec_m_y_mean = db_m.opt_fit_r_y.mean;
        s.T_rec_m_o_mean = db_m.opt_fit_r_o.mean;
        s.T_rec_y_mean = db_t.opt_fit_r_y.mean;
        s.T_rec_o_mean = db_t.opt_fit_r_o.mean;
        s.time_r = reshape(db_s.opt_fit_r_y.time_grid(1:s.k_rec+1),1,[]);
        s.k_sick = 25;
        s.time_k = reshape((0:s.k_sick),1,[]);
        s.obj_ir_y = makedist('Gamma','a',s.T_sick_y*s.T_sick_std^2,'b',1/s.T_sick_std^2);
        s.pdf_ir_y = pdf(s.obj_ir_y,s.time_k);
        s.pdf_ir_y = s.pdf_ir_y./sum(s.pdf_ir_y);
        s.T_rec_i_y_mean = dot(s.time_k,s.pdf_ir_y);
        s.obj_ir_o = makedist('Gamma','a',s.T_sick_o*s.T_sick_std^2,'b',1/s.T_sick_std^2);
        s.pdf_ir_o = pdf(s.obj_ir_o,s.time_k);
        s.pdf_ir_o = s.pdf_ir_o./sum(s.pdf_ir_o);
        s.T_rec_i_o_mean = dot(s.time_k,s.pdf_ir_o);
        s.obj_sick_y = makedist('Gamma','a',s.T_sick_y*s.T_sick_std^2,'b',1/s.T_sick_std^2);
        s.obj_sick_o = makedist('Gamma','a',s.T_sick_o*s.T_sick_std^2,'b',1/s.T_sick_std^2);
    end

    function [y]=cut_tail(y,k)
        n=length(y);
        y(end-k:end-1) = NaN;
        y(end)=0;
        y = interp1(find(~isnan(y)),y(find(~isnan(y))),1:n,'spline'); %#ok<FNDSB>
        y = y'/sum(y);
    end
end

