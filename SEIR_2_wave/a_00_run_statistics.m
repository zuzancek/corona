%% cleanup
initialize;
s = struct;

%% load data
filename = 'data/calibration/hosp_data.xlsx';
tb = readtable(filename);
save('data/calibration/hosp_data_raw.mat','tb');

%% initialization
x = tb.day;
N_rand = 50000;
s_width = 10; s_type = 5; s_ends = 1;
% s_pre = 5;

%% calculations
% 1. Hospital admission (time to, probability)
s.alpha_h_y = tb.P_H_Y(end);
s.alpha_h_o = tb.P_H_O(end);
cdf_h_y = tb.P_H_Y./s.alpha_h_y;
pdf_h_y = cdf_h_y(2:end)-cdf_h_y(1:end-1);  pdf_h_y = [pdf_h_y;0]/sum(pdf_h_y);
pdf_h_y_s = smooth_series(pdf_h_y);pdf_h_y_s = pdf_h_y_s/sum(pdf_h_y_s);
cdf_h_o = tb.P_H_O./s.alpha_h_o;
pdf_h_o = cdf_h_o(2:end)-cdf_h_o(1:end-1);  pdf_h_o = [pdf_h_o;0]/sum(pdf_h_o);
pdf_h_o_s = smooth_series(pdf_h_o);pdf_h_o_s = pdf_h_o_s/sum(pdf_h_o_s);
s.opt_fit_h_y = get_prob_dist(N_rand,x,pdf_h_y_s,'do_plot',true,'title','PHY');
s.opt_fit_h_o = get_prob_dist(N_rand,x,pdf_h_o_s,'do_plot',true,'title','PHO');
% 2. Death (time to, probability ... during hospital stay)
s.alpha_d_y = tb.P_D_Y(end);
s.alpha_d_o = tb.P_D_O(end);  
[pdf_d_y_s,cdf_d_y_s] = apply_censoring(tb.P_D_Y./s.alpha_d_y,25);
[pdf_d_o_s,cdf_d_o_s] = apply_censoring(tb.P_D_O./s.alpha_d_o,25);
s.opt_fit_d_y = get_prob_dist(N_rand,x,pdf_d_y_s,'do_plot',true,'title','PDY');
s.opt_fit_d_o = get_prob_dist(N_rand,x,pdf_d_o_s,'do_plot',true,'title','PDO');
% 3. Recovery (time to, probability ... during hospital stay)
s.alpha_r_y = tb.P_R_Y(end);
s.alpha_r_o = tb.P_R_O(end);  
[pdf_r_y_s,cdf_r_y_s] = apply_censoring(tb.P_R_Y./s.alpha_r_y,25);
[pdf_r_o_s,cdf_r_o_s] = apply_censoring(tb.P_R_O./s.alpha_r_o,25);
s.opt_fit_r_y = get_prob_dist(N_rand,x,pdf_r_y_s,'do_plot',true,'title','PRY');
s.opt_fit_r_o = get_prob_dist(N_rand,x,pdf_r_o_s,'do_plot',true,'title','PRO');

%% saving stuff
save('results/optimal_fit.mat','s');
