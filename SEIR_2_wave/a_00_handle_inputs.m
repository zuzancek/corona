initialize;

x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility.csv','dateFormat','yyyy-mm-dd','freq','daily');
hosp = dbload('data/hospitals.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_age = dbload('data/new_cases_age.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_deaths = dbload('data/deaths.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_deaths_age = dbload('data/age_deaths_cases.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_asympt = dbload('data/asymptomatical_cases_share.csv','dateFormat','yyyy-mm-dd','freq','daily');

%% definitions
s = setparam();
disp_from = dd(2020,9,1);
cut = 0;
dt = 1;
t0 = startdate(x.ActiveCases);
tt0 = t0+dt;
t1 = enddate(x.ActiveCases)-cut;
dateFrom = dd(2020,9,1);
dateTo = t1;
out_filename = 'results/inputs_full.mat';

%% data processing
% epidemiological
y = process_inputs(x,tt0,t1);
[cases_data,test_data] = process_epidemiology_data(x,y,tt0,t1,s);

% clinical
[hosp_data,deaths_data] = process_clinical_statistics(hosp,db_deaths,db_deaths_age,dateFrom-s.firstData_offset-s.k_death-1,dateTo);

%% calculations
% asymptomatic share
asymp_ratio = db_asympt.Net;
[asymp_ratio,asymp_ratio_smooth,asymp_ratio_raw] = extend_series(asymp_ratio,t0,t1,[],[]);
cases_data.asymp_ratio = asymp_ratio;      cases_data.asymp_ratio_smooth = asymp_ratio_smooth; cases_data.asymp_ratio_raw = asymp_ratio_raw;

% old-age share (in cases, dead)
old_ratio = db_age.Old./db_age.Total;
[old_ratio,old_ratio_smooth,old_ratio_raw] = extend_series(old_ratio,t0,t1,s.old_share,[]);
cases_data.old_ratio = old_ratio;          cases_data.old_ratio_smooth = old_ratio_smooth;     cases_data.old_ratio_raw = old_ratio_raw;

%% plotting
% clinical statistics
plot_clinical_statistics(hosp_data+deaths_data,dateFrom,dateTo,'raw',false,'smooth',true,'mm',true);
% mobility
mob_data = plot_mobility(mob,dateFrom,dateTo);
% epidemiology
plot_epidemiology_reporting(cases_data,test_data,dateFrom,dateTo,'raw',false,'mm',true,'smooth',true);

%% saving
dates.t0 = t0;      dates.t1 = disp_from;   dates.t2 = t1;
save(out_filename,'dates','cases_data','test_data','hosp_data','deaths_data','mob_data','s');