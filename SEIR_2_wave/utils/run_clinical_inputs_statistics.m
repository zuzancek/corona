function [s]=run_clinical_inputs_statistics(filename)

tb = readtable(filename);
save('data/calibration/hosp_data_raw.mat','tb');

%% initialization
x = tb.day;
N_rand = 50000;
s_width = 10; s_type = 5; s_ends = 1;
last_idx_d = 25;
last_idx_r = 25;

%% calculations
% 1. Hospital admission (time to, probability)
alpha_h_y = tb.P_H_Y(end);
alpha_h_o = tb.P_H_O(end);
cdf_h_y = tb.P_H_Y./alpha_h_y;
pdf_h_y = cdf_h_y(2:end)-cdf_h_y(1:end-1);  pdf_h_y = [pdf_h_y;0]/sum(pdf_h_y);
pdf_h_y_s = smooth_series(pdf_h_y);pdf_h_y_s = pdf_h_y_s/sum(pdf_h_y_s);
cdf_h_o = tb.P_H_O./alpha_h_o;
pdf_h_o = cdf_h_o(2:end)-cdf_h_o(1:end-1);  pdf_h_o = [pdf_h_o;0]/sum(pdf_h_o);
pdf_h_o_s = smooth_series(pdf_h_o);pdf_h_o_s = pdf_h_o_s/sum(pdf_h_o_s);
s.opt_fit_h_y = get_prob_dist(N_rand,x,pdf_h_y_s,'do_plot',true,'title','PHY');
s.opt_fit_h_y.alpha = alpha_h_y;
s.opt_fit_h_o = get_prob_dist(N_rand,x,pdf_h_o_s,'do_plot',true,'title','PHO');
s.opt_fit_h_o.alpha = alpha_h_o;
% 2. Death (time to, probability ... during hospital stay)
alpha_d_y = tb.P_D_Y(end);
alpha_d_o = tb.P_D_O(end);  
[pdf_d_y_s,~] = apply_censoring(tb.P_D_Y./alpha_d_y,last_idx_d,s_width,s_type,s_ends);
[pdf_d_o_s,~] = apply_censoring(tb.P_D_O./alpha_d_o,last_idx_d,s_width,s_type,s_ends);
s.opt_fit_d_y = get_prob_dist(N_rand,x,pdf_d_y_s,'do_plot',true,'title','PDY');
s.opt_fit_d_y.alpha = alpha_d_y;
s.opt_fit_d_o = get_prob_dist(N_rand,x,pdf_d_o_s,'do_plot',true,'title','PDO');
s.opt_fit_d_o.alpha = alpha_d_o;
% 3. Recovery (time to, probability ... during hospital stay)
alpha_r_y = tb.P_R_Y(end);
alpha_r_o = tb.P_R_O(end);  
[pdf_r_y_s,~] = apply_censoring(tb.P_R_Y./alpha_r_y,last_idx_r,s_width,s_type,s_ends);
[pdf_r_o_s,~] = apply_censoring(tb.P_R_O./alpha_r_o,last_idx_r,s_width,s_type,s_ends);
s.opt_fit_r_y = get_prob_dist(N_rand,x,pdf_r_y_s,'do_plot',true,'title','PRY');
s.opt_fit_r_y.alpha = alpha_r_y;
s.opt_fit_r_o = get_prob_dist(N_rand,x,pdf_r_o_s,'do_plot',true,'title','PRO');
s.opt_fit_r_o.alpha = alpha_r_o;

end
