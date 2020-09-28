initialize;

x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
s = setparam();

%% handle data
dt = 1;
t0 = startdate(x.ActiveCases);
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
% q_mat_new = zeros(length(s.quant), 6+length(Rt));
% q_mat_new(:,7:end) = q_mat; 
% q_mat = q_mat_new;
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
xticks([0 30 60 90 120 150 180 210])
xticklabels({'marec','apr�l','m�j','j�n','j�l','august','september'})
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
