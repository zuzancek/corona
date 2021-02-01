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

%% calculations
% hospital admission
s.alpha_h_y = tb.P_H_Y(end);
s.alpha_h_o = tb.P_H_O(end);
cdf_h_y = tb.P_H_Y./s.alpha_h_y;
pdf_h_y = cdf_h_y(2:end)-cdf_h_y(1:end-1);  pdf_h_y = [pdf_h_y;0]/sum(pdf_h_y);
pdf_h_y_s = smooth_series(pdf_h_y);pdf_h_y_s = pdf_h_y_s/sum(pdf_h_y_s);
cdf_h_o = tb.P_H_O./s.alpha_h_o;
pdf_h_o = cdf_h_o(2:end)-cdf_h_o(1:end-1);  pdf_h_o = [pdf_h_o;0]/sum(pdf_h_o);
pdf_h_o_s = smooth_series(pdf_h_o);pdf_h_o_s = pdf_h_o_s/sum(pdf_h_o_s);
s.T_h_o = dot(x,pdf_h_o)+1;                 s.T_h_y = dot(x,pdf_h_y)+1;
s.T_h_o_s = dot(x,pdf_h_o_s)+1;             s.T_h_y_s = dot(x,pdf_h_y_s)+1;
opt_fit_h_y = get_prob_dist(N_rand,x,pdf_h_y_s,'do_plot',true);
opt_fit_h_o = get_prob_dist(N_rand,x,pdf_h_o_s,'do_plot',true);

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