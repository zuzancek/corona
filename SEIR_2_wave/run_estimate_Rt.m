initialize;

x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility_new.csv','dateFormat','yyyy-mm-dd','freq','daily');

s = setparam();
disp_from = dd(2020,4,1);
indiff = true; 

%% handle data
% start date: 6.3
cut = 0;
dt = 1;
t0 = startdate(x.ActiveCases);
disp_to = enddate(x.ActiveCases)-cut;
tt0 = t0+dt;
t1 = enddate(x.ActiveCases)-cut;

dI_inflow = resize(x.NewCases,tt0:t1);
dI_inflow_smooth = smooth_series(double(dI_inflow),s.smooth_width,...
    s.smooth_type,s.smooth_ends);

dI_inflow_smooth = 0*dI_inflow+dI_inflow_smooth;

pos_test_ratio = x.NewCases./x.Tests;
pos_test_ratio_smooth = smooth_series(double(pos_test_ratio),s.smooth_width,...
    s.smooth_type,s.smooth_ends);
[dI_inflow_adj_smooth,dI_inflow_adj] = adjust_series(double(dI_inflow),s.smooth_width,...
    s.smooth_type,s.smooth_ends,s.scale_fact,s.threshold,double(resize(pos_test_ratio,tt0:t1)));
dI_inflow_adj = 0*dI_inflow+dI_inflow_adj;
dI_inflow_adj_smooth = 0*dI_inflow+dI_inflow_adj_smooth;
I0 = x.TotalCases(tt0-1)/s.obs_ratio;

%% calculations
[Rt,~,~,Xt] = estimate_Rt(double(dI_inflow),I0,s.pop_size,s.T_rem,s.sim_num);
[Rt_smooth,q_mat,It_smooth,Xt_smooth,x_mat,Rt_last,St_smooth] = estimate_Rt(double(dI_inflow_smooth),I0,s.pop_size,s.T_rem,s.sim_num,s.quant,s.pweight);
% forecast_Rt(Rt_smooth,double(dI_inflow_smooth), 50,I0,s.pop_size,s.T_rem,s.sim_num,s.quant);

[Rt_adj,~,~,Xt_adj] = estimate_Rt(double(dI_inflow_adj),I0,s.pop_size,s.T_rem,s.sim_num);
[Rt_adj_smooth,q_mat_adj,~,Xt_adj_smooth,x_mat_adj] = estimate_Rt(double(dI_inflow_adj_smooth),I0,s.pop_size,s.T_rem,s.sim_num,s.quant);

pos_test_ratio_smooth = 0*pos_test_ratio+pos_test_ratio_smooth;

%% plotting stuff
threshold = 100;
figure;
fn = fieldnames(mob);
for i=1:length(fn)
    yy = interp(resize(mob.(fn{i}),t0:t1),t0:t1);
    if i==length(fn)
        plot(yy,'linewidth',2,'Color','k');hold on;
    else
        plot(yy,'linewidth',1);hold on;
    end
end
xl = xlim;
bench = tseries(xl(1):xl(2),0)+threshold;
plot(bench,'color',[0.15 0.15 0.15],'linestyle','--');
grid on;
ylabel('%');
title('Mobility, raw data');
legend(fn);

figure;
for i=1:length(fn)
    yy = interp(resize(mob.(fn{i}),t0:t1),t0:t1);
    zz = smooth_series(double(yy),s.smooth_width,s.smooth_type,s.smooth_ends);
    yy = 0*yy+zz;
    if i==length(fn)
        plot(yy,'linewidth',2,'Color','k');hold on;
    else
        plot(yy,'linewidth',1);hold on;
    end
end
plot(bench,'color',[0.15 0.15 0.15],'linestyle','--');
grid on;
ylabel('%');
title('Mobility, smooth data');
legend(fn);

figure;
subplot(2,1,1);
plot(dI_inflow,'linewidth',1);hold on;
plot(dI_inflow_smooth,'linewidth',2);hold on;
title('New infections (observed)');
legend({'raw','smooth'});
grid on;

subplot(2,1,2);
plot(diff(dI_inflow_smooth),'linewidth',2);hold on;
plot(diff(dI_inflow_adj_smooth),'linewidth',2,'linestyle','--');hold on;
title('Change in infections (smooth)');
legend({'observed','adjusted '});
grid on;

figure;
plot(dI_inflow,'linewidth',1);hold on;
plot(dI_inflow_smooth,'linewidth',2);hold on;
plot(dI_inflow_adj,'linewidth',1,'linestyle','--');hold on;
plot(dI_inflow_adj_smooth,'linewidth',2,'linestyle','--');hold on;
title('New infections (adjusted)');
legend({'observed, raw','observed, smooth', 'hypothetical, raw','hypothetical,smooth'});
grid on;

%
figure;
plot(pos_test_ratio,'linewidth',1); hold on;
plot(pos_test_ratio_smooth,'linewidth',1);
title('Positive tests ratio');
legend({'raw','smooth'});
grid on;
%
figure;
subplot(2,1,1);
plot(0*dI_inflow+Rt,'linewidth',1);hold on;
plot(0*dI_inflow+Rt_smooth,'linewidth',1);hold on;
title('Rt');
legend({'raw','smooth'});
grid on;
subplot(2,1,2);
plot(0*dI_inflow+Rt_adj,'linewidth',1);hold on;
plot(0*dI_inflow+Rt_adj_smooth,'linewidth',1);hold on;
title('Rt (with assumption)');
legend({'raw','smooth'});
grid on;
figure;
plot(0*dI_inflow+Rt_smooth,'linewidth',1);hold on;
plot(0*dI_inflow+Rt_adj_smooth,'linewidth',1);hold on;
title('Rt smooth');
legend({'observed','adjusted'});
grid on;
% 
plot_fanchart(q_mat,s,dt,disp_from,disp_to,t0);
plot_fanchart(q_mat_adj,s,dt,disp_from,disp_to,t0);
plot_fanchart(x_mat,s,dt,disp_from,disp_to,t0);
plot_fanchart(x_mat_adj,s,dt,disp_from,disp_to,t0);

%% saving stuff
Rt_vec_raw = zeros(t1-t0+1,1);
Rt_vec_raw(dt+1:end) = Rt;
Rt_vec_raw = tseries(t0:t1,Rt_vec_raw);
x.Rt_raw = Rt_vec_raw;

Rt_vec_smooth = zeros(t1-t0+1,1);
Rt_vec_smooth(dt+1:end) = Rt_smooth;
Rt_vec_smooth = tseries(t0:t1,Rt_vec_smooth);
x.Rt_smooth = Rt_vec_smooth;

dbsave(x,'results.csv');
dIt = dI_inflow_smooth;
St = St_smooth;
It = It_smooth;
Rt = Rt_smooth;
save('inputs.mat','q_mat','Rt','dIt','It','St','s','Rt_last','t0','t1');
