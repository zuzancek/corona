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
dateFrom = dd(2020,9,1);
dateTo = t1;

% epidemiological
y = process_inputs(x,tt0,t1);
[cases_data,test_data] = process_epidemiology_data(x,y,tt0,t1,s);

% clinical
[hosp_data,deaths_data] = process_clinical_statistics(hosp,db_deaths,db_deaths_age,dateFrom-s.firstData_offset,dateTo);

%% calculations
% asymptomatic share
asymp_ratio = db_asympt.Net;
[asymp_ratio,asymp_ratio_smooth,asymp_ratio_raw] = extend_series(asymp_ratio,t0,t1,[],[]);
cases_data.asymp_ratio = asymp_ratio;      cases_data.asymp_ratio_smooth = asymp_ratio_smooth; cases_data.asymp_ratio_raw = asymp_ratio_raw;

% old-age share (in cases, dead)
old_ratio = db_age.Old./db_age.Total;
[old_ratio,old_ratio_smooth,old_ratio_raw] = extend_series(old_ratio,t0,t1,s.old_share,[]);
cases_data.old_ratio = old_ratio;          cases_data.old_ratio_smooth = old_ratio_smooth;     cases_data.old_ratio_raw = old_ratio_raw;

% observed ratio
delay.v = 1*([1 2 1 0]);  delay.at = [dd(2020,9,15),dd(2020,10,15),dd(2020,11,15),dd(2020,12,15)];
params = struct;
params.death_old_ratio = deaths_data.old_ratio;
params.cases_old_ratio = old_ratio;
params.asymp_ratio = asymp_ratio;
params.cutoff = 3;
params.adj = 0; % *0
params.h = hosp_data.H;
[dI_inflow_real, I_real, obs_ratio_real,sa_cmp,par] = fun_0(x,hosp_data,s,disp_from,t1,t0,t1,...
    params,delay);
cases_data.cases_pcr_implied = dI_inflow_real;
cases_data.cases_pcr_implied_smooth = smooth_series(dI_inflow_real);
cases_data.obs_ratio = obs_ratio_real;
cases_data.loss = sa_cmp;
cases_data.X_smooth = par.X_smooth;
cases_data.X_forecast = par.X_forecast_smooth;
cases_data.X_total = tseries(startdate(par.X_smooth):enddate(par.X_forecast_smooth),0);
cases_data.X_total(startdate(par.X_smooth):enddate(par.X_smooth)) = par.X_smooth;
cases_data.X_total(enddate(par.X_smooth)+1:enddate(par.X_forecast_smooth)) = par.X_forecast_smooth;
cases_data.cases_pcr_implied_smooth = cases_data.X_total;

% alternative numbers for hospitals
init = hosp_data;
init.I = x.ActiveCases; 
init.rho = old_ratio; init.varsigma = db_deaths_age.TotalDeathRatioOld;
[out] = fun_1(x,par,s,init,disp_from,t1);
y = x; y.NewCases = dI_inflow_real;
[out_check] = fun_1(y,par,s,init,disp_from,t1);
hosp_data.alt = out;

%% plotting stuff
% 1. reporting
% clinical statistics
plot_clinical_statistics(hosp_data+deaths_data,dateFrom,dateTo,'raw',false,'smooth',true,'mm',true);
% mobility
mob_data = plot_mobility(mob,dateFrom,dateTo);
% epidemiology
plot_epidemiology_reporting(cases_data,test_data,dateFrom,dateTo,'raw',false,'mm',true,'smooth',true);

% 2. analysis, comparison
% epidemiology
% hospitals
plot_clinical_cmp(out_check,hosp_data+deaths_data,out,dateFrom,dateTo,'mm',true,'reported',true,'reduced',true);

figure('Name','New cases (reported vs.true)');
fh1 = plot(par.X_rep_smooth,'linewidth',2);hold on;
fh2 = plot(par.X_smooth,'linewidth',3);hold on; 
fh3 = plot(resize(cases_data.cases_total_smooth,disp_from:t1),'linewidth',2);
% fh3 = plot(resize(dI_inflow_smooth,disp_from:t1),'linewidth',2);hold on;
plot(par.X_forecast_smooth,'linewidth',3, 'linestyle',':','Color',fh2.Color);
% plot(par.X_raw,'Color',[0.75 0.75 0.75],'linewidth',1); 
% plot(par.X_forecast_raw,'Color',[0.55 0.55 0.55],'linewidth',1,'linestyle',':');
% plot(par.X_rep_raw,'linewidth',1,'Color',[0.5 0.5 0.5],'linestyle','-.');
% plot(resize(cases_data.cases_pcr, startdate(par.X_rep_raw):enddate(par.X_rep_raw)),'linewidth',1,'Color',[0.75 0.75 0.75],'linestyle','-.');
% plot(resize(dI_inflow,disp_from:t1),'linewidth',1,'Color',[0.5 0.5 0.5]);
d_from = startdate(par.X_rep_raw);
d_to = enddate(par.X_rep_raw);
lab = {'september','oktober','november','december','januar'};
mtt = d_from+cumsum([0 30 31 30 31]);
xticks(mtt)
xticklabels(lab);
ylim([0 8000]);
xlim([d_from dd(2021,1,29)]);
legend([fh1 fh2 ],{'Reported (confirmed) new cases (PCR tests)','Implied by hospitals/deaths (+forecast)', 'Reported cases (PCR+AG)'});%, 'Reported new cases (PCR+AG tests)'}); 
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
xls_out.test_eff = 100/s.obs_ratio*resize(obs_ratio_real,disp_from:t1);
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

%% saving stuff
dates.t0 = t0;      dates.t1 = disp_from;   dates.t2 = t1;
save(out_filename,'dates','cases_data','test_data','hosp_data','deaths_data','mob_data','s');