initialize;

x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility.csv','dateFormat','yyyy-mm-dd','freq','daily');
hosp = dbload('data/hospitals.csv','dateFormat','yyyy-mm-dd','freq','daily');

s = setparam();
disp_from = dd(2020,4,1);
indiff = true; 

%% handle data
% definitions
cut = 0;
dt = 1;
t0 = startdate(x.ActiveCases);
disp_to = enddate(x.ActiveCases)-cut;
tt0 = t0+dt;
t1 = enddate(x.ActiveCases)-cut;
% t1 = dd(2020,11,14);
h_t0 = startdate(hosp.ICU);
h_t00 = find(hosp.ICU>0);
h_t1 = enddate(hosp.ICU);

% epidemiology data
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

% clinical
hospit = hosp.Hospitalizations;
hospit_smooth = 0*hospit+smooth_series(double(hospit),5,...
    s.smooth_type,s.smooth_ends);
icu = hosp.ICU;
icu_smooth = 0*resize(icu,h_t00:h_t1)+smooth_series(icu(h_t00:h_t1),5,...
    s.smooth_type,s.smooth_ends);
vent = hosp.Ventilation;
vent_smooth = 0*vent+smooth_series(double(vent),5,...
    s.smooth_type,s.smooth_ends);
death = resize(x.Deaths,h_t0:h_t1);
death_smooth = 0*death+smooth_series(double(death),5,...
    s.smooth_type,s.smooth_ends);

%% calculations
pos_test_ratio_smooth = 0*pos_test_ratio+pos_test_ratio_smooth;

%% plotting stuff
% clinical statistics
figure('Name','Clinical statistics');
subplot(2,2,1)
plot(hospit,'linestyle','--');hold on;
plot(hospit_smooth,'linewidth',1);
grid on;
title('Total hospitalizations');
subplot(2,2,2)
plot(icu,'linestyle','--');hold on;
plot(icu_smooth,'linewidth',1);
grid on;
title('ICU beds');
subplot(2,2,3)
plot(vent,'linestyle','--');hold on;
plot(vent_smooth,'linewidth',1);
grid on;
title('Ventilations');
subplot(2,2,4)
plot(death,'linestyle','--');hold on;
plot(death_smooth,'linewidth',1);
grid on;
title('Deaths');

% mobility
threshold = 100;
figure('Name','Mobility');
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

% epidemiology
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

figure;
plot(pos_test_ratio,'linewidth',1); hold on;
plot(pos_test_ratio_smooth,'linewidth',1);
title('Positive tests ratio');
legend({'raw','smooth'});
grid on;

%% saving stuff
mob_smooth = yy;
save('inputs.mat','dI_inflow','dI_inflow_smooth','dI_inflow_adj','dI_inflow_adj_smooth',...
    'pos_test_ratio','pos_test_ratio_smooth','I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1','h_t00');
data.NewCases = x.NewCases;
data.Tests = x.Tests;
data.Mobility = mob.SK;
data.Deaths = x.Deaths;
data.ICU = icu;
data.Vent = vent;
data.Hosp = hospit;
save('raw.mat','data');