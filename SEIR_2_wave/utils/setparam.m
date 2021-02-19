function [s] = setparam(varargin)

s = struct;
s.pop_size = 5443120;
s.sim_num = 10^4;
s.T = 100;
s.model_seir = true;
s.tests_min = 4700;
s.cases_min = 750;
s.ratio_threshold = 0.05;
s.env_from = dd(2020,10,8);
s.wave_2_from = dd(2020,9,1);
s.firstData_offset = 31;
s.deaths_with_covid_share = .2;
s.S_H_rate = 0.16;
s.S_H_rate_0 = 0.25;
s.scale_s_h = s.S_H_rate_0/s.S_H_rate;

% serial interval (generation period)
s.SI.mean = 7.5;                s.SI.std = 0.62;
% time to test (observation period, from symptoms onset): "steady_state value"
s.T_test.mean = 2;              s.T_test.std = s.SI.std;      
% incubation period 
s.T_inc.mean = 5.2;             s.T_inc.std = s.SI.std;
% infectious period
s.T_inf.mean = 4.3;             s.T_inf.std = 0.62;
s.T_inf_asymp.mean = 4.3;       s.T_inf_asymp.std = 0.62;
s.T_inf_symp.mean = 4.3;        s.T_inf_symp.std = 0.62;
s.T_inf_obs.mean = 4.3;         s.T_inf_obs.std = 0.62;
s.T_inf_unobs.mean = 4.3;       s.T_inf_unobs.std = 0.62;
% sickness/symptoms period
s.T_sick_y = 11;                s.T_sick_o = 15;       s.T_sick = 13;
s.T_sick_std = s.SI.std;
s.k_sick = 20;                  s.T_sick_pdf_type = 'Gamma'; 
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

try
    db = load('results/optimal_fit.mat','stat_total','stat_severe','stat_mild');
    set_prob_data();
catch err %#ok<NASGU>
    a_00_run_statistics();
    db = load('results/optimal_fit.mat','s');
    set_prob_data();
end

% total time shift in clinical model
s.t_shift_clin = 30;

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

s.smooth_width = 7;
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

s.smoothing_method_data = @mov_median_adj;
s.smoothing_method_params = @smooth_series;

    function[]=set_prob_data()
        % dbs definitions
        db_t = db.stat_total; db_s = db.stat_severe; db_m = db.stat_mild;
        % death (serious cases only, use db_s database)
        s.k_death = 40;
        s.omega_y = db_s.opt_fit_d_y.alpha; %/(1-s.deaths_with_covid_share);
        s.omega_o = db_s.opt_fit_d_o.alpha; %/(1-s.deaths_with_covid_share);
        s.pdf_d_y = max(0,db_s.opt_fit_d_y.pdf(1:s.k_death+1))/sum(max(0,db_s.opt_fit_d_y.pdf(1:s.k_death+1)));
        s.pdf_d_o = max(0,db_s.opt_fit_d_o.pdf(1:s.k_death+1))/sum(max(0,db_s.opt_fit_d_o.pdf(1:s.k_death+1)));
        s.T_death_y_mean = db_s.opt_fit_d_y.mean;
        s.T_death_o_mean = db_s.opt_fit_d_o.mean;
        s.time_d = db_s.opt_fit_d_y.time_grid(1:s.k_death+1);
        % admission to ICU, ventilation, ECMO
        s.k_ser = 15;
        s.theta_y = db_s.opt_fit_h_y.alpha;
        s.theta_o = db_s.opt_fit_h_o.alpha;
        s.pdf_s_y = db_s.opt_fit_h_y.pdf(1:s.k_ser+1)/sum(db_s.opt_fit_h_y.pdf(1:s.k_ser+1));
        s.pdf_s_o = db_s.opt_fit_h_o.pdf(1:s.k_ser+1)/sum(db_s.opt_fit_h_o.pdf(1:s.k_ser+1));
        s.T_ser_y_mean = db_s.opt_fit_h_y.mean;
        s.T_ser_o_mean = db_s.opt_fit_h_o.mean;
        % hospital admission
        s.k_hosp = 20;        
        s.eta_y = db_t.opt_fit_h_y.alpha*1.5;
        s.eta_o = db_t.opt_fit_h_o.alpha*0.75;        
        s.pdf_m_y = db_t.opt_fit_h_y.pdf(1:s.k_hosp+1)/sum(db_t.opt_fit_h_y.pdf(1:s.k_hosp+1));
        s.pdf_m_o = db_t.opt_fit_h_o.pdf(1:s.k_hosp+1)/sum(db_t.opt_fit_h_o.pdf(1:s.k_hosp+1));        
        s.T_hosp_y_mean = db_t.opt_fit_h_y.mean;
        s.T_hosp_o_mean = db_t.opt_fit_h_o.mean;
        s.time_h = db_t.opt_fit_h_y.time_grid(1:s.k_hosp+1);
        % recovery (serious,moderate,home)
        s.k_rec = 40;
        s.pdf_sr_y = db_s.opt_fit_r_y.pdf(1:s.k_rec+1)/sum(db_s.opt_fit_r_y.pdf(1:s.k_rec+1));
        s.pdf_sr_o = db_s.opt_fit_r_o.pdf(1:s.k_rec+1)/sum(db_s.opt_fit_r_o.pdf(1:s.k_rec+1));
        s.T_rec_s_y_mean = db_s.opt_fit_r_y.mean;
        s.T_rec_s_o_mean = db_s.opt_fit_r_o.mean;
        s.pdf_mr_y = db_m.opt_fit_r_y.pdf(1:s.k_rec+1)/sum(db_m.opt_fit_r_y.pdf(1:s.k_rec+1));
        s.pdf_mr_o = db_m.opt_fit_r_o.pdf(1:s.k_rec+1)/sum(db_m.opt_fit_r_o.pdf(1:s.k_rec+1));
        s.T_rec_m_y_mean = db_m.opt_fit_r_y.mean;
        s.T_rec_m_o_mean = db_m.opt_fit_r_o.mean;
        s.T_rec_y_mean = db_t.opt_fit_r_y.mean;
        s.T_rec_o_mean = db_t.opt_fit_r_o.mean;
        s.time_r = db_s.opt_fit_r_y.time_grid(1:s.k_rec+1);
        s.k_sick = 20;
        s.time_s = (0:s.k_sick)';
        s.pdf_ir_y = pdf('Gamma',s.time_s,s.T_sick_y*s.T_sick_std^2,1/s.T_sick_std^2); 
        s.pdf_ir_y = s.pdf_ir_y./sum(s.pdf_ir_y);        
        s.T_rec_i_y_mean = dot(s.time_s,s.pdf_ir_y);
        s.pdf_ir_o = pdf('Gamma',s.time_s,s.T_sick_o*s.T_sick_std^2,1/s.T_sick_std^2); 
        s.pdf_ir_o = s.pdf_ir_o./sum(s.pdf_ir_o);
        s.T_rec_i_o_mean = dot(s.time_s,s.pdf_ir_o); 
        s.obj_sick_y = makedist('Gamma','a',s.T_sick_y*s.T_sick_std^2,'b',1/s.T_sick_std^2);
        s.obj_sick_o = makedist('Gamma','a',s.T_sick_o*s.T_sick_std^2,'b',1/s.T_sick_std^2);        
    end

%     function[]=set_prob_data()
%         % hospital admission
%         s.k_hosp = 20;
%         db_t = db.stat_total; db_s = db.stat_severe; db_m = db.stat_mild;
%         s.eta_y = db_t.opt_fit_h_y.alpha*1.25;
%         s.pdf_h_y = db_t.opt_fit_h_y.pdf(1:s.k_hosp+1)/sum(db_t.opt_fit_h_y.pdf(1:s.k_hosp+1));
%         s.T_hosp_y_mean = db_t.opt_fit_h_y.mean;
%         s.time_h = db_t.opt_fit_h_y.time_grid(1:s.k_hosp+1);
%         s.eta_o = db_t.opt_fit_h_o.alpha*0.75;
%         s.pdf_h_o = db_t.opt_fit_h_o.pdf(1:s.k_hosp+1)/sum(db_t.opt_fit_h_o.pdf(1:s.k_hosp+1));
%         s.T_hosp_o_mean = db_t.opt_fit_h_o.mean;
%         s.obj_hosp_y = db_t.opt_fit_h_y.obj;
%         s.obj_hosp_o = db_t.opt_fit_h_o.obj;
%         % death
%         s.k_death = 40;
%         s.omega_y = db_t.opt_fit_d_y.alpha/(1-s.deaths_with_covid_share); 
%         s.omega_y_m = db_m.opt_fit_d_y.alpha/(1-s.deaths_with_covid_share); 
%         s.omega_y_s = db_s.opt_fit_d_y.alpha/(1-s.deaths_with_covid_share);
%         s.pdf_d_y = db_t.opt_fit_d_y.pdf(1:s.k_death+1)/sum(db_t.opt_fit_d_y.pdf(1:s.k_death+1));
%         s.T_death_o_mean = db_t.opt_fit_d_o.mean;
%         s.time_d = db_t.opt_fit_d_y.time_grid(1:s.k_death+1);
%         s.omega_o = db_t.opt_fit_d_o.alpha/(1-s.deaths_with_covid_share); 
%         s.omega_o_m = db_m.opt_fit_d_o.alpha/(1-s.deaths_with_covid_share); 
%         s.omega_o_s = db_s.opt_fit_d_o.alpha/(1-s.deaths_with_covid_share); 
%         s.pdf_d_o = db_t.opt_fit_d_o.pdf(1:s.k_death+1)/sum(db_t.opt_fit_d_o.pdf(1:s.k_death+1));
%         s.T_death_y_mean = db_t.opt_fit_d_y.mean;
%         % recovery at hospital
%         s.k_rec = 40;
%         s.pdf_r_y = db_t.opt_fit_r_y.pdf(1:s.k_rec+1)/sum(db_t.opt_fit_r_y.pdf(1:s.k_rec+1));      
%         s.T_rec_y_mean = db_t.opt_fit_r_y.mean;
%         s.T_rec_y_s_mean = db_s.opt_fit_r_y.mean;
%         s.pdf_r_y_s = db_s.opt_fit_r_y.pdf(1:s.k_rec+1)/sum(db_s.opt_fit_r_y.pdf(1:s.k_rec+1));  
%         s.obj_r_y_s = db_s.opt_fit_r_y.obj;
%         s.T_rec_y_m_mean = db_m.opt_fit_r_y.mean;
%         s.pdf_r_y_m = db_m.opt_fit_r_y.pdf(1:s.k_rec+1)/sum(db_m.opt_fit_r_y.pdf(1:s.k_rec+1)); 
%         s.obj_r_y = db_t.opt_fit_r_y.obj;
%         s.T_rec_o_mean = db_t.opt_fit_r_o.mean;
%         s.time_r = db_t.opt_fit_r_y.time_grid(1:s.k_rec+1);
%         s.pdf_r_o = db_t.opt_fit_r_o.pdf(1:s.k_rec+1)/sum(db_t.opt_fit_r_o.pdf(1:s.k_rec+1));
%         s.T_rec_o_s_mean = db_s.opt_fit_r_o.mean;
%         s.pdf_r_o_s = db_s.opt_fit_r_o.pdf(1:s.k_rec+1)/sum(db_s.opt_fit_r_o.pdf(1:s.k_rec+1));
%         s.obj_r_o_s = db_s.opt_fit_r_o.obj;
%         s.T_rec_o_m_mean = db_m.opt_fit_r_o.mean;
%         s.pdf_r_o_m = db_m.opt_fit_r_o.pdf(1:s.k_rec+1)/sum(db_m.opt_fit_r_o.pdf(1:s.k_rec+1));
%         s.obj_r_o = db_t.opt_fit_r_o.obj;
%         % recovery at home
%         s.k_sick = 20;
%         s.time_s = (0:s.k_sick)';
%         s.pdf_s_y = pdf('Gamma',s.time_s,s.T_sick_y*s.T_sick_std^2,1/s.T_sick_std^2); 
%         s.pdf_s_y = s.pdf_s_y./sum(s.pdf_s_y);        
%         s.T_sick_y_mean = dot(s.time_s,s.pdf_s_y);
%         s.pdf_s_o = pdf('Gamma',s.time_s,s.T_sick_o*s.T_sick_std^2,1/s.T_sick_std^2); 
%         s.pdf_s_o = s.pdf_s_o./sum(s.pdf_s_o);
%         s.T_sick_o_mean = dot(s.time_s,s.pdf_s_o); 
%         s.obj_sick_y = makedist('Gamma','a',s.T_sick_y*s.T_sick_std^2,'b',1/s.T_sick_std^2);
%         s.obj_sick_o = makedist('Gamma','a',s.T_sick_o*s.T_sick_std^2,'b',1/s.T_sick_std^2);        
%     end
end

