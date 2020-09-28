initialize;

x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
s = setparam();

%% handle data
% first case: 6.3.2020
% T0 = 18.3.2020  ... 12 days
% T1 = 1.6.2020 ..... 88 days
dt = 1;
t0 = startdate(x.ActiveCases);
%t1 = t0+88-12-1;
tt0 = t0+dt;
t1 = enddate(x.ActiveCases);
dI_inflow = x.NewCases(tt0:t1);
dI_inflow_smooth = smooth_series(dI_inflow,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
I0 = x.TotalCases(tt0-1)/s.obs_ratio;
[Rt] = estimate_Rt(dI_inflow,I0,s.pop_size,s.T_rem,s.sim_num);
[Rt_smooth,q_mat] = estimate_Rt(dI_inflow_smooth,I0,s.pop_size,s.T_rem,s.sim_num,s.quant);

%% plotting stuff
figure;
plot(dI_inflow,'linewidth',1);hold on;
plot(dI_inflow_smooth,'linewidth',1);hold on;
title('Wave 1: New infections, observed');
legend({'raw','smooth'});
grid on;

figure;
plot(Rt,'linewidth',1);hold on;
plot(Rt_smooth,'linewidth',1);hold on;
title('Wave 1: Rt');
legend({'raw','smooth'});
grid on;

f = figure;
q_mat = q_mat(:,1:end);
fanChart(1:size(q_mat,2), q_mat', q_mat(s.quant_idx_central,:), s.quant,...
    'alpha', .75, 'colormap', {'shadesOfColor',s.color_graph});
grid on;
set(gca,'color',s.color_bkg);
% axis tight
ax = gca;
ax.GridColor = s.color_grid;
ax.FontName = 'TimesNewRoman';
ax.FontWeight = 'bold';
ax.XAxis.Color = s.color_grid;
ax.YAxis.Color = s.color_grid;
f.Color = s.color_bkg;
% legend(s.quant_legend);


%% savig stuff
Rt_vec_raw = zeros(t1-t0+1,1);
Rt_vec_raw(dt+1:end) = Rt;
Rt_vec_raw = tseries(t0:t1,Rt_vec_raw);
x.Rt_raw = Rt_vec_raw;

Rt_vec_smooth = zeros(t1-t0+1,1);
Rt_vec_smooth(dt+1:end) = Rt_smooth;
Rt_vec_smooth = tseries(t0:t1,Rt_vec_smooth);
x.Rt_smooth = Rt_vec_smooth;

dbsave(x,'results.csv');
