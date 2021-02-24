initialize;

x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility.csv','dateFormat','yyyy-mm-dd','freq','daily');
hosp = dbload('data/hospitals.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_deaths = dbload('data/deaths.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_deaths_age = dbload('data/age_deaths_cases.csv','dateFormat','yyyy-mm-dd','freq','daily');

load('results/inputs.mat','dates','cases_data','test_data','hosp_data','deaths_data','s');

idx_fun = 1;
fun_opt_0 = {'DHIX','DHIXt'}; fun_0 = str2func(fun_opt_0{idx_fun});
fun_opt_1 = {'XIHD','XIHDt'}; fun_1 = str2func(fun_opt_1{idx_fun});
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

%% calculations
% 1./ daily cases implied by hospitals
delay.v = 1*([1 0.5 0]);  delay.at = [dd(2020,10,30),dd(2020,11,15),dd(2020,12,15)];
params = struct;
params.death_old_ratio = deaths_data.old_ratio;
params.cases_old_ratio = cases_data.old_ratio;
params.asymp_ratio = cases_data.asymp_ratio;
params.death_adj = deaths_data.delta;
params.serious_cases_ratio = hosp_data.S_H_rate;
params.cutoff = 3;
params.adj = 0; 
params.h = hosp_data.H;
params.s = hosp_data.S;
[dI_inflow_real, res_implied, params] = fun_0(x,hosp_data,s,disp_from,t1,t0,params,delay);
cases_data.cases_pcr_implied = dI_inflow_real;
cases_data.cases_pcr_implied_smooth = smooth_series(dI_inflow_real);
cases_data.obs_ratio = res_implied.obs_ratio_adj;
cases_data.loss = res_implied.sa;
cases_data.X_smooth = res_implied.X_smooth;
cases_data.X_forecast = res_implied.X_forecast_smooth;
cases_data.X_total = res_implied.X_smooth_total;
cases_data.cases_pcr_implied_smooth = cases_data.X_total;

% 2./ hospitals/deaths implied by reported cases 
init = hosp_data;
init.I = x.ActiveCases; 
init.T_delay = params.T_delay;
init.IH = res_implied.IH;
init.HR = res_implied.HR;
init.omega_o = s.omega_o; 
init.omega_y = s.omega_y;
init.rho = old_ratio; 
init.varsigma = db_deaths_age.TotalDeathRatioOld;
[out] = fun_1(x,params,s,init,disp_from,t1);

% 3./ hospitals/deaths implied by implied cases 
% backcheck technique
y = x; 
y.NewCases = res_implied.X_all;
init.rho = params.rho;
init.kappa_d = params.kappa_d;
init.nu = params.nu;
init.kappa_h_o = params.kappa_h_o;
init.kappa_h_y = params.kappa_h_y;
init.kappa_s = params.kappa_s;
init.omega_o = params.omega_o; init.omega_y = params.omega_y;
[out_check] = fun_1(y,params,s,init,disp_from,t1);
hosp_data.alt = out;

%% plotting stuff
% hospitals
reported_data = hosp_data+deaths_data;
true_data = out_check;
counterfact_data = out;
plot_clinical_cmp(true_data,counterfact_data,reported_data,dateFrom,dateTo,'mm',true,'reported',true);

% epidemiology
plot_epidemiology_cmp(res_implied,counterfact_data,reported_data,dateFrom,dateTo,'implied',true,'reported',true);


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