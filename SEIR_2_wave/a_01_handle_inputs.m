initialize;

x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility.csv','dateFormat','yyyy-mm-dd','freq','daily');
hosp = dbload('data/hospitals.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_age = dbload('data/new_cases_age.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_deaths = dbload('data/deaths.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_deaths_age = dbload('data/age_deaths_cases.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_asympt = dbload('data/asymptomatical_cases_share.csv','dateFormat','yyyy-mm-dd','freq','daily');

s = setparam();
idx_fun = 2;
out_filename_opt = {'results/inputs_full.mat','results/inputs.mat'}; out_filename = out_filename_opt{idx_fun};
fun_opt_0 = {'DHIXe','DHIXt'}; fun_0 = str2func(fun_opt_0{idx_fun});
fun_opt_1 = {'XIHDe','XIHDt'}; fun_1 = str2func(fun_opt_1{idx_fun});
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
dI_inflow_pcr(end) = dI_inflow_pcr(end-4);
dI_inflow_pcr_smooth = smooth_series(s.smoothing_method_data(dI_inflow_pcr));
dI_inflow = dI_inflow_pcr+dI_inflow_ag;
dI_inflow_smooth = smooth_series(s.smoothing_method_data(dI_inflow));

pos_test_ratio = x.NewCases./x.Tests;
pos_test_ratio_ag = y.AgPosit./y.AgTests;
pos_test_ratio_smooth = smooth_series(s.smoothing_method_data(pos_test_ratio));
tests = x.Tests;
tests_smooth = smooth_series(s.smoothing_method_data(tests));
I0 = x.TotalCases(t0)/s.obs_ratio;

cases_data = struct;
cases_data.cases_ag = dI_inflow_ag;   cases_data.cases_ag_mm = mov_median(dI_inflow_ag);   cases_data.cases_ag_smooth = smooth_series(dI_inflow_ag);
cases_data.cases_pcr = dI_inflow_pcr; cases_data.cases_pcr_mm = mov_median(dI_inflow_pcr); cases_data.cases_pcr_smooth = dI_inflow_pcr_smooth;
cases_data.cases_total = dI_inflow;   cases_data.cases_total_mm = mov_median(dI_inflow);   cases_data.cases_total_smooth = dI_inflow_smooth;
cases_data.ptr_pcr = pos_test_ratio;  cases_data.ptr_mm = mov_median(pos_test_ratio);      cases_data.ptr_smooth = pos_test_ratio_smooth;
cases_data.ptr_ag = pos_test_ratio_ag;cases_data.ptr_ag_mm = mov_median(pos_test_ratio_ag);cases_data.ptr_ag_smooth = tests_smooth;
cases_data.I0 = I0;

% clinical
hospit = hosp.Hospitalizations;
hospit_smooth = smooth_series(hospit);
icu = hosp.ICU;
icu(startdate(hospit_smooth)+40:startdate(hospit_smooth)+60) = NaN;
icu = interp(icu,startdate(hospit_smooth):enddate(hospit_smooth));
hosp.ICU = icu;
icu_smooth = smooth_series(resize(icu,h_t0:h_t1),s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
vent = hosp.Ventilation;
vent_smooth = smooth_series(vent,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
death = resize(x.Deaths,h_t0:h_t1);
death_smooth = smooth_series(death,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
admiss = hosp.Admission;
admiss_smooth = smooth_series(admiss,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
discharge = hosp.Discharge;
discharge_smooth = smooth_series(discharge,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
deaths_total = db_deaths.Total;
deaths_total_smooth = smooth_series(deaths_total,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
deaths_onCovid = db_deaths.DeathCovid;
deaths_onCovid_smooth = smooth_series(deaths_onCovid,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);
deaths_withCovid = db_deaths.DeathWithCovid;
deaths_withCovid_smooth = smooth_series(deaths_withCovid,s.smooth_width_hosp,s.smooth_type,s.smooth_ends);

hosp_data = struct;
hosp_data.H_raw = hospit;           hosp_data.H_smooth = hospit_smooth;
hosp_data.C_raw = icu;              hosp_data.C_smooth = icu_smooth;
hosp_data.V_raw = vent;             hosp_data.V_smooth = vent_smooth;
hosp_data.D_raw = deaths_total;     hosp_data.D_smooth = deaths_total_smooth;

deaths_data = struct;
deaths_data.total = deaths_total; deaths_data.total_smooth = deaths_total_smooth;
deaths_data.on = deaths_onCovid; deaths_data.on_smooth = deaths_onCovid_smooth;
deaths_data.with = deaths_withCovid; deaths_data.with_smooth = deaths_withCovid_smooth;

%% calculations
% asymptomatic share
asymp_ratio = db_asympt.Net;
[asymp_ratio,asymp_ratio_smooth] = extend_series(asymp_ratio,t0,t1,[],[]);
cases_data.asymp_ratio = asymp_ratio;      cases_data.asymp_ratio_smooth = asymp_ratio_smooth;

% old-age share (in cases, dead)
old_ratio = db_age.Old./db_age.Total;
[old_ratio,old_ratio_smooth] = extend_series(old_ratio,t0,t1,s.old_share,[]);
cases_data.old_ratio = old_ratio;          cases_data.old_ratio_smooth = old_ratio_smooth;

% case fatality rate at hospitals
cfr_init = []; cfr_final = 17.5;
[cfr,cfr_smooth,cfr_ext,cfr_ext_smooth] = process_xls('data/cfr_hospitals.xlsx', dd(2020,10,15),dd(2020,12,08),dd(2020,3,13),t1,s,cfr_init,cfr_final);
deaths_data.cfr = cfr_ext;                  deaths_data.cfr_smooth = cfr_ext_smooth;

% observed ratio
delay.v = [1 2 0.5 0];  delay.at = [dd(2020,10,1),dd(2020,10,31),dd(2020,11,15),dd(2020,12,15)];
% delay.v = [0.25 0.75 0];  delay.at = [dd(2020,10,25),dd(2020,11,15),dd(2020,12,15)];
params = struct;
params.death_old_ratio = db_deaths_age.TotalDeathRatioOld;
deaths_data.old_ratio = params.death_old_ratio; 
deaths_data.old_ratio_smooth = smooth_series(deaths_data.old_ratio);
params.cfr_hospitals = cfr_ext;
params.cases_old_ratio = old_ratio;
params.asymp_ratio = asymp_ratio;
other.mob = mob.SK;
other.ptr = pos_test_ratio;
params.other = other;
params.cutoff = 4;
params.firstData = disp_from-31;
params.adj = -0.175;
[dI_inflow_real, I_real, obs_ratio_real,sa_cmp,par] = fun_0(x,hosp,deaths_total,s,disp_from,t1,t0,t1,...
    params,delay);
cases_data.cases_pcr_implied = dI_inflow_real;
cases_data.cases_pcr_implied_smooth = smooth_series(dI_inflow_real);
cases_data.obs_ratio = obs_ratio_real;
cases_data.loss = sa_cmp;

% alternative numbers for hospitals
init.D = death; init.H = hospit; 
init.C = icu;   init.V = vent;
init.I = x.ActiveCases; 
init.rho = old_ratio; init.varsigma = db_deaths_age.TotalDeathRatioOld;
[out] = fun_1(x,par,s,init,disp_from,t1);
y = x; y.NewCases = dI_inflow_real;
[out_check] = fun_1(y,par,s,init,disp_from,t1);
hosp_data.alt = out;

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
pp1=plot(resize(death,disp_from:t1),'linewidth',1,'linestyle','--');hold on;
pp2=plot(resize(deaths_total,disp_from:t1),'linewidth',1,'linestyle','--');hold on;
plot(resize(death_smooth,disp_from:t1),'linewidth',2,'Color',pp1.Color);
plot(resize(deaths_total_smooth,disp_from:t1),'linewidth',2,'Color',pp2.Color);
grid on;
title('Total Deaths (with+on covid)');

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

figure('Name','Mobility, smooth data');
for i=1:length(fn)
    yy = interp(resize(mob.(fn{i}),t0:t1),t0:t1);
    zz = smooth_series(yy,s.smooth_width,s.smooth_type,s.smooth_ends);
    if i==length(fn)
        plot(zz,'linewidth',2,'Color','k');hold on;
    else
        plot(zz,'linewidth',1);hold on;
    end
end
plot(bench,'color',[0.15 0.15 0.15],'linestyle','--');
grid on;
ylabel('%');
title('Mobility, smooth data');
legend(fn);

% epidemiology
figure('Name','New cases, PCR, AG')
plot(resize(dI_inflow_pcr,disp_from:t1),'linewidth',1);hold on;
plot(resize(dI_inflow_pcr_smooth,disp_from:t1),'linewidth',2);hold on;
plot(resize(dI_inflow,disp_from:t1),'linewidth',1);hold on;
plot(resize(dI_inflow_smooth,disp_from:t1),'linewidth',2);hold on;
title('New infections (reported only)');
legend({'PCR only: raw','PCR only: smooth','PCR+AG: raw','PCR+AG: smooth'});
grid on;

figure('Name','CFR & asymptomatic share')
subplot(2,1,1);
plot(resize(100*cfr,disp_from:t1),'linewidth',1);hold on;
plot(resize(100*cfr_ext_smooth,disp_from:t1),'linewidth',2);hold on;
title('Case fatality rate (reported,  %, hospitals only)');
legend({'raw','smooth'});
grid on;

subplot(2,1,2);
plot(resize(100*asymp_ratio,disp_from:t1),'linewidth',1);hold on;
plot(resize(100*asymp_ratio_smooth,disp_from:t1),'linewidth',2);hold on;
title('Share of asymptomatic new cases (reported,  %, PCR only)');
legend({'raw','smooth'});
grid on;

figure('Name','Tests & PTR');
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

figure('Name','New cases (reported vs.true)');
fh1 = plot(par.X_rep_smooth,'linewidth',2);hold on;
fh2 = plot(par.X_smooth,'linewidth',3);hold on; 
% fh3 = plot(resize(dI_inflow_smooth,disp_from:t1),'linewidth',2);hold on;
plot(par.X_forecast_smooth,'linewidth',3, 'linestyle',':','Color',fh2.Color);
plot(par.X_raw,'Color',[0.55 0.55 0.55],'linewidth',1); 
% plot(par.X_forecast_raw,'Color',[0.55 0.55 0.55],'linewidth',1,'linestyle',':');
plot(par.X_rep_raw,'linewidth',1,'Color',[0.5 0.5 0.5]);
% plot(resize(dI_inflow,disp_from:t1),'linewidth',1,'Color',[0.5 0.5 0.5]);
grid on;
legend([fh1 fh2 ],{'Reported (confirmed) new cases (PCR tests)','Implied by hospitals/deaths (+forecast)'});%, 'Reported new cases (PCR+AG tests)'}); 
title('New cases (smooth data)');

figure('Name','New cases (reported vs.true, lost cases)');
plot(resize(dI_inflow_pcr,disp_from:t1),'linewidth',1,'linestyle','-.');hold on;
plot(resize(dI_inflow_pcr_smooth,disp_from:t1),'linewidth',2);hold on;
plot(resize(sa_cmp.loss_a,disp_from:t1),'linewidth',1);hold on;
plot(resize(dI_inflow_real,disp_from:t1),'linewidth',2);hold on;
plot(resize(sa_cmp.loss_s,disp_from:t1),'linewidth',1);hold on;
title('New infections (PCR only)');
legend({'reported, raw','reported, smooth', '"lost" asymptomatical new cases','hypothetically observable','"lost" symptomatical new cases'});
grid on;
% 
figure('Name','Testing effectivity and Old-age cases share')
subplot(2,1,1)
obs_ratio = 0*obs_ratio_real+s.obs_ratio;
plot(100*resize(obs_ratio,disp_from:t1),'linewidth',1); hold on;
plot(100*resize(obs_ratio_real,disp_from:t1),'linewidth',1);
title('Observable ratio');
legend({'stationary (optimistic)','real (implied by hospitalizations)'});
grid on;
ylabel('% of total cases');
subplot(2,1,2)
plot(100/s.obs_ratio*resize(obs_ratio_real,disp_from:t1),'linewidth',1);grid on;
title('Testing effectivity (implied by hospitals)');

% 
figure('Name','Testing effectivity and Old-age cases share')
subplot(2,1,1)
obs_ratio = 0*obs_ratio_real+s.obs_ratio;
plot(100*resize(obs_ratio,disp_from:t1),'linewidth',1); hold on;
plot(100*resize(obs_ratio_real,disp_from:t1),'linewidth',1);
title('Observable ratio');
legend({'stationary (optimistic)','real (implied by hospitalizations)'});
grid on;
ylabel('% of total cases');
subplot(2,1,2)
plot(resize(old_ratio,disp_from:t1),'linewidth',1); hold on;
plot(resize(old_ratio_smooth,disp_from:t1),'linewidth',1); hold on;
title('Old-age persons (65+, % of confirmend cases)');
legend({'raw','smooth'});
grid on;
ylabel('%');

%
figure('Name','Situation in Hospitals: Comparison I')
if idx_fun==1
    subplot(2,2,1)
    plot(resize(hospit_smooth,disp_from:t1),'linewidth',1);hold on;
    plot(resize(out.H,disp_from:t1),'linewidth',1);hold on;
    legend({'observed','implied by reported daily new cases'});
    grid on;
    title('Hospitalisations (total)');
    subplot(2,2,2)
    plot(resize(icu_smooth,disp_from:t1),'linewidth',1);hold on;
    plot(resize(out.C,disp_from:t1),'linewidth',1);hold on;
    grid on;
    title('ICU');
    subplot(2,2,3)
    plot(resize(vent_smooth,disp_from:t1),'linewidth',1);hold on;
    plot(resize(out.V,disp_from:t1),'linewidth',1);hold on;
    grid on;
    title('Ventilations');
    subplot(2,2,4)
    plot(resize(deaths_total_smooth,disp_from:t1),'linewidth',1);hold on;
    plot(resize(out.D,disp_from:t1),'linewidth',1);hold on;
    grid on;
    title('Deaths');
else
    subplot(2,1,1)
    ratio_h = resize(hospit_smooth,disp_from:t1)./resize(out_check.H,disp_from:t1);
    plot(resize(hospit_smooth,disp_from:t1),'linewidth',2);hold on;
    plot(ratio_h.*resize(out.H,disp_from:t1),'linewidth',2);hold on;
    % plot(ratio.*resize(out_check.H,disp_from:t1),'k--','linewidth',1);hold on;
    legend({'observed','implied by reported daily new cases'});%,'reconstructed from implied cases'});
    grid on;
    title('Hospitalisations (total)');
    subplot(2,1,2)
    ratio_d = resize(deaths_total_smooth,disp_from:t1)./resize(out_check.D,disp_from:t1);
    plot(resize(deaths_total_smooth,disp_from:t1),'linewidth',2);hold on;
    plot(ratio_d.*resize(out.D,disp_from:t1),'linewidth',2);hold on;
    % plot(resize(out_check.D,disp_from:t1),'k--','linewidth',1);hold on;
    grid on;
    title('Deaths');
    legend({'observed','implied by reported daily new cases'}); %,'reconstructed from implied cases'}); 
    figure('Name','Situation in Hospitals: Comparison II');
    subplot(2,1,1)
    ratio_c = 1/3*ratio_d+2/3*ratio_h;
    plot(resize(icu_smooth,disp_from:t1),'linewidth',2);hold on;
    plot(resize(icu_smooth,disp_from:t1)./(2/3*resize(hospit_smooth,disp_from:t1)./(ratio_h.*resize(out.H,disp_from:t1))+1/3*resize(deaths_total_smooth,disp_from:t1)./(ratio_d.*resize(out.D,disp_from:t1))),'linewidth',2);hold on;
    legend({'observed','implied by reported daily new cases'});%,'reconstructed from implied cases'});
    grid on;
    title('ICU');
    subplot(2,1,2)
    plot(resize(vent_smooth,disp_from:t1),'linewidth',2);hold on;
    plot(resize(vent_smooth,disp_from:t1)./(1/3*resize(hospit_smooth,disp_from:t1)./(ratio_h.*resize(out.H,disp_from:t1))+2/3*resize(deaths_total_smooth,disp_from:t1)./(ratio_d.*resize(out.D,disp_from:t1))),'linewidth',2);hold on;
    grid on;
    title('Ventilations');
    legend({'observed','implied by reported daily new cases'}); %,'reconstructed from implied cases'}); 
end    

%% saving stuff
dates.t0 = t0;      dates.t1 = disp_from;   dates.t2 = t1;
mob_data.raw = yy;  mob_data.smooth = zz;
save(out_filename,'dates','cases_data','hosp_data','deaths_data','mob_data','s');
