% initialize;

% y = dbload('data/antigen_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility.csv','dateFormat','yyyy-mm-dd','freq','daily');
hosp = dbload('data/hospitals.csv','dateFormat','yyyy-mm-dd','freq','daily');

s = setparam();
disp_from = dd(2020,9,1);
indiff = true; 

%% handle data
% definitions
cut = 0;
dt = 1;
t0 = startdate(x.ActiveCases);
t0_ag = startdate(x.AgTests);
disp_to = enddate(x.ActiveCases)-cut;
tt0 = t0+dt;
t1 = enddate(x.ActiveCases)-cut-0;
t1_ag = enddate(x.AgTests)-cut-0;
h_t0 = startdate(hosp.ICU);
h_t00 = find(hosp.ICU>0);
h_t1 = enddate(hosp.ICU);

% epidemiology data
y = process_inputs(x,tt0,t1);
dI_inflow_ag = y.AgPosit;
dI_inflow_pcr = resize(x.NewCases,tt0:t1);
dI_inflow_pcr_smooth = smooth_series(dI_inflow_pcr,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
dI_inflow = dI_inflow_pcr+dI_inflow_ag;
dI_inflow_smooth = smooth_series(dI_inflow,s.smooth_width,...
    s.smooth_type,s.smooth_ends);

pos_test_ratio = x.NewCases./x.Tests;
pos_test_ratio_ag = y.AgPosit./y.AgTests;
pos_test_ratio_smooth = smooth_series(pos_test_ratio,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
tests = x.Tests;
tests_smooth = smooth_series(tests,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
I0 = x.TotalCases(tt0-1)/s.obs_ratio;

% clinical
hospit = hosp.Hospitalizations;
hospit_smooth = smooth_series(hospit,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
icu = hosp.ICU;
icu_smooth = smooth_series(resize(icu,h_t00:h_t1),s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
vent = hosp.Ventilation;
vent_smooth = smooth_series(vent,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
death = resize(x.Deaths,h_t0:h_t1);
death_smooth = smooth_series(death,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);

%% calculations
% observed ratio
[obs_ratio_smooth,obs_ratio,dI_inflow_pcr_adj_smooth,dI_inflow_pcr_adj,t_tests] = adjust_observed_ratio(...
    pos_test_ratio_smooth,dI_inflow_pcr_smooth,s,t0);

% asymptomatic share
final.date = t1; final.value = 15;
initial.date = t0; initial.value = 25; 
try
    [asymp_ratio,asymp_ratio_smooth] = process_as('data/asympt_share.xlsx',dd(2020,3,13),dd(2020,10,13),s,initial,final);
catch err
    pr = load('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_pcr_adj','dI_inflow_smooth','dI_inflow_pcr_adj_smooth',...
        'pos_test_ratio_smooth','obs_ratio_smooth','asymp_ratio_smooth',...
        'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
        'death_smooth','h_t0','h_t1','h_t00');
    asymp_ratio_smooth = pr.asymp_ratio_smooth;
    pr = load('raw.mat','data');
    asymp_ratio = pr.data.Asymp;
end
    
%% plotting stuff
% clinical statistics
figure('Name','Clinical statistics');
subplot(2,2,1)
plot(hospit,'linewidth',1);hold on;
plot(hospit_smooth,'linewidth',2);
grid on;
title('Total hospitalizations');
subplot(2,2,2)
plot(icu,'linewidth',1);hold on;
plot(icu_smooth,'linewidth',2);
grid on;
title('ICU beds');
subplot(2,2,3)
plot(vent,'linewidth',1);hold on;
plot(vent_smooth,'linewidth',2);
grid on;
title('Ventilations');
subplot(2,2,4)
plot(death,'linewidth',1);hold on;
plot(death_smooth,'linewidth',2);
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
plot(resize(dI_inflow_pcr,disp_from:t1),'linewidth',1);hold on;
plot(resize(dI_inflow_pcr_smooth,disp_from:t1),'linewidth',2);hold on;
plot(resize(dI_inflow,disp_from:t1),'linewidth',1);hold on;
plot(resize(dI_inflow_smooth,disp_from:t1),'linewidth',2);hold on;
title('New infections (observed)');
legend({'PCR only: raw','PCR only: smooth','PCR+AG: raw','PCR+AG: smooth'});
grid on;

subplot(2,1,2);
plot(resize(100*asymp_ratio,disp_from:t1),'linewidth',1);hold on;
plot(resize(100*asymp_ratio_smooth,disp_from:t1),'linewidth',2);hold on;
title('Share of asymptomatic new cases (observed,  %, PCR only)');
legend({'raw','smooth'});
grid on;

figure;
subplot(2,1,1)
plot(resize(pos_test_ratio,disp_from:t1),'linewidth',1); hold on;
plot(resize(pos_test_ratio_smooth,disp_from:t1),'linewidth',2);
title('Positive tests ratio');
legend({'PCR: raw','PCR: smooth'});
grid on;
subplot(2,1,2);
plot(resize(tests,disp_from:t1),'linewidth',1);hold on;
plot(resize(tests_smooth,disp_from:t1),'linewidth',2);hold on;
title('Tests (observed)');
legend({'PCR: raw','PCR: smooth'});
grid on;

figure;
plot(resize(dI_inflow_pcr,disp_from:t1),'linewidth',1,'linestyle','-.');hold on;
plot(resize(dI_inflow_pcr_smooth,disp_from:t1),'linewidth',2);hold on;
plot(resize(dI_inflow_pcr_adj,disp_from:t1),'linewidth',1,'linestyle','-.');hold on;
plot(resize(dI_inflow_pcr_adj_smooth,disp_from:t1),'linewidth',2);hold on;
plot(resize(dI_inflow_smooth,disp_from:t1),'linewidth',2);hold on;
title('New infections (PCR only)');
legend({'observed, raw','observed, smooth', 'hypothetical, raw','hypothetical,smooth','AG included in observed'});
grid on;

figure;
subplot(2,1,1)
plot(resize(obs_ratio,disp_from:t1),'linewidth',1); hold on;
plot(resize(obs_ratio_smooth,disp_from:t1),'linewidth',1);
title('Hypothetical observed ratio');
legend({'PCR: raw','PCR: smooth'});
grid on;

%% saving stuff
mob_smooth = yy;
obs_ratio_smooth = resize(obs_ratio_smooth,startdate(dI_inflow_pcr_smooth):enddate(dI_inflow_pcr_smooth));
asymp_ratio_smooth = resize(asymp_ratio_smooth,startdate(dI_inflow_pcr_smooth):enddate(dI_inflow_pcr_smooth));
pos_test_ratio_smooth = resize(pos_test_ratio_smooth,startdate(dI_inflow_pcr_smooth):enddate(dI_inflow_pcr_smooth));
save('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_pcr_adj','dI_inflow_smooth','dI_inflow_pcr_adj_smooth',...
    'pos_test_ratio_smooth','obs_ratio_smooth','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','h_t0','h_t1','h_t00');
data.NewCases_PCR = x.NewCases;
data.Tests_PCR = x.Tests;
data.Asymp = asymp_ratio;
data.Obs_PCR = obs_ratio;
data.TestRatio_PCR = pos_test_ratio;
data.Mobility = mob.SK;
data.Deaths = x.Deaths;
data.ICU = icu;
data.Vent = vent;
data.Hosp = hospit;
save('raw.mat','data');