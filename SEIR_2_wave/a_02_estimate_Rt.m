%% initialization
initialize;

%% setup
disp_from = dd(2020,4,1);
indiff = true; 
cut = 0;
dt = 1;


%% load data
load('inputs.mat','dI_inflow','dI_inflow_smooth','dI_inflow_adj','dI_inflow_adj_smooth',...
    'pos_test_ratio','pos_test_ratio_smooth','I0','mob','s','t0','t1');
tt0 = t0+dt;
t1 = dd(2020,11,14);

if s.model_seir
    model_fnc = @estimate_Rt_SEIR_aug;
    del = 1;
else
    model_fnc = @estimate_Rt_SIR;
    del = 0;
end
disp_to = t1-del-1;

%% calculations
[Rt,~,~,Xt] = model_fnc(double(resize(dI_inflow_smooth,t0:t1)),I0,s);
[Rt_smooth,q_mat,It_smooth,Xt_smooth,x_mat,Rt_last,St_smooth,Et_smooth] = model_fnc(double(resize(dI_inflow_smooth,t0:t1)),I0,s,s.quant,s.pweight);

[Rt_adj,~,~,Xt_adj] = model_fnc(double(resize(dI_inflow_adj,t0:t1)),I0,s);
[Rt_adj_smooth,q_mat_adj,~,Xt_adj_smooth,x_mat_adj] = model_fnc(double(resize(dI_inflow_adj_smooth,t0:t1)),I0,s,s.quant);

pos_test_ratio_smooth = 0*pos_test_ratio+pos_test_ratio_smooth;

%% plotting stuff
y = 0*resize(dI_inflow,t0:t1-del);
figure;
subplot(2,1,1);
plot(y+Rt,'linewidth',1);hold on;
plot(y+Rt_smooth,'linewidth',1);hold on;
title('Rt');
legend({'raw','smooth'});
grid on;
subplot(2,1,2);
plot(y+Rt_adj,'linewidth',1);hold on;
plot(y+Rt_adj_smooth,'linewidth',1);hold on;
title('Rt (with assumption)');
legend({'raw','smooth'});
grid on;
figure;
plot(y+Rt_smooth,'linewidth',1);hold on;
plot(y+Rt_adj_smooth,'linewidth',1);hold on;
title('Rt smooth');
legend({'observed','adjusted'});
grid on;
% 
plot_fanchart(q_mat,s,dt,disp_from,disp_to,t0,'Effective reproduction number (Rt)');
plot_fanchart(q_mat_adj,s,dt,disp_from,disp_to,t0,'Adjusted effective reproduction number (Rt)');
plot_fanchart(x_mat,s,dt,disp_from,disp_to,t0,'Active infections estimate (unobs.included)');
plot_fanchart(x_mat_adj,s,dt,disp_from,disp_to,t0,'Adjusted active infections estimate (unobs.included)');

%% saving stuff
Rt_vec_raw = zeros(t1-t0+1-del,1);
Rt_vec_raw(dt+1:end) = Rt;
Rt_vec_raw = tseries(t0:t1-del,Rt_vec_raw);
x.Rt_raw = Rt_vec_raw;

Rt_vec_smooth = zeros(t1-t0+1-del,1);
Rt_vec_smooth(dt+1:end) = Rt_smooth;
Rt_vec_smooth = tseries(t0:t1-del,Rt_vec_smooth);
x.Rt_smooth = Rt_vec_smooth;

dbsave(x,'results.csv');
dIt = dI_inflow_smooth;
St = tseries(t0+1:t1,St_smooth);
It = tseries(t0+1:t1,It_smooth);
Et = tseries(t0+1:t1,Et_smooth);
Rt = tseries(t0+1:t1-del,Rt_smooth);
save('results_Rt.mat','q_mat','Rt','dIt','It','Et','St','s','Rt_last','t0','t1');
