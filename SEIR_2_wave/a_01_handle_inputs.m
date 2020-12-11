initialize;

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
icu_smooth = smooth_series(resize(icu,h_t0:h_t1),s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
vent = hosp.Ventilation;
vent_smooth = smooth_series(vent,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
death = resize(x.Deaths,h_t0:h_t1);
death_smooth = smooth_series(death,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
admiss = hosp.Admission;
admiss_smooth = smooth_series(admiss,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
discharge = hosp.Discharge;
discharge_smooth = smooth_series(discharge,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);

%% calculations
% asymptomatic share
final.date = t1; final.value = 13.5;
initial.date = t0; initial.value = 25; 
breakpoint.date = dd(2020,10,30);breakpoint.value = 29.5;
try
    [asymp_ratio,asymp_ratio_smooth] = process_as('data/asympt_share.xlsx',dd(2020,3,13),dd(2020,10,30),...
        s,initial,final,breakpoint,dd(2020,11,30));
catch err
    pr = load('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_real','dI_inflow_smooth','dI_inflow_real_smooth',...
        'pos_test_ratio_smooth','obs_ratio','asymp_ratio_smooth',...
        'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
        'death_smooth','h_t0','h_t1');
    asymp_ratio_smooth = pr.asymp_ratio_smooth;
    pr = load('raw.mat','data');
    asymp_ratio = pr.data.Asymp;
end

% observed ratio
[dI_inflow_real, I_real, obs_ratio_real,sa_cmp] = adjust_infection_hospitals_full(x,hosp,s,disp_from,t1,t0,t1,asymp_ratio_smooth);
    
%% plotting stuff
% clinical statistics
figure('Name','Clinical statistics');
subplot(2,2,1)
plot(resize(hospit,disp_from:t1),'linewidth',1);hold on;
plot(resize(hospit_smooth,disp_from:t1),'linewidth',2);
grid on;
title('Total hospitalizations');
subplot(2,2,2)
plot(resize(icu,disp_from:t1),'linewidth',1);hold on;
plot(resize(icu_smooth,disp_from:t1),'linewidth',2);
grid on;
title('ICU beds');
subplot(2,2,3)
plot(resize(vent,disp_from:t1),'linewidth',1);hold on;
plot(resize(vent_smooth,disp_from:t1),'linewidth',2);
grid on;
title('Ventilations');
subplot(2,2,4)
plot(resize(death,disp_from:t1),'linewidth',1);hold on;
plot(resize(death_smooth,disp_from:t1),'linewidth',2);
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
figure('Name','New cases & asymptomatic share');
subplot(2,1,1);
plot(resize(dI_inflow_pcr,disp_from:t1),'linewidth',1);hold on;
plot(resize(dI_inflow_pcr_smooth,disp_from:t1),'linewidth',2);hold on;
plot(resize(dI_inflow,disp_from:t1),'linewidth',1);hold on;
plot(resize(dI_inflow_smooth,disp_from:t1),'linewidth',2);hold on;
title('New infections (reported only)');
legend({'PCR only: raw','PCR only: smooth','PCR+AG: raw','PCR+AG: smooth'});
grid on;

subplot(2,1,2);
plot(resize(100*asymp_ratio,disp_from:t1),'linewidth',1);hold on;
plot(resize(100*asymp_ratio_smooth,disp_from:t1),'linewidth',2);hold on;
title('Share of asymptomatic new cases (reported,  %, PCR only)');
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
title('Tests (reported)');
legend({'PCR: raw','PCR: smooth'});
grid on;

figure('Name','Observable ratio and lost cases');
subplot(2,1,1)
plot(resize(dI_inflow_pcr,disp_from:t1),'linewidth',1,'linestyle','-.');hold on;
plot(resize(dI_inflow_pcr_smooth,disp_from:t1),'linewidth',2);hold on;
plot(resize(sa_cmp.loss_a,disp_from:t1),'linewidth',1);hold on;
plot(resize(dI_inflow_real,disp_from:t1),'linewidth',2);hold on;
plot(resize(sa_cmp.loss_s,disp_from:t1),'linewidth',1);hold on;
title('New infections (PCR only)');
legend({'reported, raw','reported, smooth', '"lost" asymptomatical new cases','hypothetically observable','"lost" symptomatical new cases'});
grid on;

subplot(2,1,2)
obs_ratio = 0*obs_ratio_real+s.obs_ratio;
plot(100*resize(obs_ratio,disp_from:t1),'linewidth',1); hold on;
plot(100*resize(obs_ratio_real,disp_from:t1),'linewidth',1);
title('Observable ratio');
legend({'stationary (optimistic)','real (implied by hospitalizations)'});
grid on;
ylabel('% of total cases');


%% saving stuff
mob_smooth = yy;
% obs_ratio_smooth = resize(obs_ratio_smooth,startdate(dI_inflow_pcr_smooth):enddate(dI_inflow_pcr_smooth));
asymp_ratio_smooth = resize(asymp_ratio_smooth,startdate(dI_inflow_pcr_smooth):enddate(dI_inflow_pcr_smooth));
pos_test_ratio_smooth = resize(pos_test_ratio_smooth,startdate(dI_inflow_pcr_smooth):enddate(dI_inflow_pcr_smooth));
save('inputs.mat','dI_inflow_pcr','dI_inflow_pcr_smooth','dI_inflow_real','dI_inflow_smooth',...
    'pos_test_ratio_smooth','obs_ratio','obs_ratio_real','asymp_ratio_smooth',...
    'I0','mob','s','t0','t1','hospit_smooth','vent_smooth','icu_smooth',...
    'death_smooth','admiss_smooth','discharge_smooth','h_t0','h_t1');
data.NewCases_PCR = x.NewCases;
data.TotalActiveCases_real_PCR = I_real;
data.NewCases_real_PCR = dI_inflow_real;
data.Tests_PCR = x.Tests;
data.Asymp = asymp_ratio;
data.ObsRatio_real = obs_ratio_real;
data.TestRatio_PCR = pos_test_ratio;
data.Mobility = mob.SK;
data.Deaths = x.Deaths;
data.ICU = icu;
data.Vent = vent;
data.Hosp = hospit;
data.Admission = admiss;
data.Discharge = discharge;
save('raw.mat','data');