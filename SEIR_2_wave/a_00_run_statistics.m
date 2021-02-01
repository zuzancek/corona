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
s.opt_fit_h_y = get_prob_dist(N_rand,x,pdf_h_y_s,'competitive_risk',false,'do_plot',true,'title','PHY');
s.opt_fit_h_o = get_prob_dist(N_rand,x,pdf_h_o_s,'competitive_risk',false,'do_plot',true,'title','PHO');
% 2. Death (time to, probability ... during hospital stay)
s.alpha_d_y = tb.P_D_Y(end);
s.alpha_d_o = tb.P_D_O(end);
cdf_d_y = tb.P_D_Y./s.alpha_d_y;            
cdf_d_y_s = smooth_series(cdf_d_y,s_width,s_type,s_ends);             cdf_d_y_s(1) = 0;
pdf_d_y = cdf_d_y_s(2:end)-cdf_d_y_s(1:end-1);  pdf_d_y = [pdf_d_y;0]/sum(pdf_d_y);
pdf_d_y_s = smooth_series(pdf_d_y,s_width,s_type,s_ends);pdf_d_y_s = pdf_d_y_s/sum(pdf_d_y_s);
cdf_d_o = tb.P_D_O./s.alpha_d_o;
pdf_d_o = cdf_d_o(2:end)-cdf_d_o(1:end-1);  pdf_d_o = [pdf_d_o;0]/sum(pdf_d_o);
pdf_d_o_s = smooth_series(pdf_d_o);pdf_d_o_s = pdf_d_o_s/sum(pdf_d_o_s);
s.opt_fit_d_y = get_prob_dist(N_rand,x,pdf_d_y_s,'competitive_risk',false,'do_plot',true,'title','PDY');
s.opt_fit_d_o = get_prob_dist(N_rand,x,pdf_d_o_s,'competitive_risk',false,'do_plot',true,'title','PDO');

%% plotting 
% hospital admission
figure('Name','Hospital admission I.');
xho = find(tb.P_H_O(2:end)-tb.P_H_O(1:end-1)==0,1);
xhy = find(tb.P_H_Y(2:end)-tb.P_H_Y(1:end-1)==0,1);
xh_1 = max(xho,xhy);
xh = 0:xh_1;
subplot(2,1,1);
stairs(xh,tb.P_H_Y(1:xh_1+1),'linewidth',2);
hold on;
stairs(xh,tb.P_H_O(1:xh_1+1),'linewidth',2);
legend({'P(H(t<=T)|Y)','P(H(t<=T)|O)'});
grid on;
ylabel('%');
xlabel('days from test');
title('Endpoint/hospitalisation probability up to specific date');
%
subplot(2,1,2);
stairs(xh,100*cdf_h_y(1:xh_1+1),'linewidth',2);
hold on;
stairs(xh,100*cdf_h_o(1:xh_1+1),'linewidth',2);
legend({'P(H(t<=T)|Y)','P(H(t<=T)|O)'});
grid on;
ylabel('%');
xlabel('days from test');
title('CDF');
%
figure('Name','Hospital admission II.');
xho = find(pdf_h_o(2:end)-pdf_h_o(1:end-1)==0,1);
xhy = find(pdf_h_y(2:end)-pdf_h_y(1:end-1)==0,1);
xh_1 = max(xho,xhy)-1;
xh = 0:xh_1;
subplot(2,1,1);
stairs(xh,pdf_h_y(1:xh_1+1),'-.','linewidth',1);
hold on;
plot(xh,pdf_h_y_s(1:xh_1+1),'linewidth',2);
grid on;
ylabel('%');
xlabel('days from test');
title('Prob(H|Y)');

subplot(2,1,2);
stairs(xh,pdf_h_o(1:xh_1+1),'-.','linewidth',1);
hold on;
plot(xh,pdf_h_o_s(1:xh_1+1),'linewidth',2);
grid on;
ylabel('%');
xlabel('days from test');
title('Prob(H|O)');