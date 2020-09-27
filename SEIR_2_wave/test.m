initialize;

x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
s = setparam();

%% First wave
% first case: 6.3.2020
% T0 = 18.3.2020  ... 12 days
% T1 = 1.6.2020 ..... 88 days
t0 = startdate(x.ActiveCases)+12;
t1 = t0+88-12-1;
dI_inflow = x.NewCases(t0:t1);
dI_inflow_smooth = smooth_series(dI_inflow,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
I0 = x.TotalCases(t0-1)/s.obs_ratio;
[Rt,Rt_vec] = estimate_Rt(dI_inflow,I0,s.pop_size,s.T_rem);
[Rt_smooth,Rt_vec_smooth] = estimate_Rt(dI_inflow_smooth,I0,s.pop_size,s.T_rem);

figure;
plot(dI_inflow,'linewidth',1);hold on;
plot(dI_inflow_smooth,'linewidth',1);hold on;
title('Wave 1: New infections, observed');
legend({'raw','smooth'});
grid on;

