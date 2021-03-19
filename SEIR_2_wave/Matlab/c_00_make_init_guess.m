
%% initialization
initialize;

%% load data
db = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_age = dbload('data/new_cases_age.csv','dateFormat','yyyy-mm-dd','freq','daily');
s = setparam();

% 1 = reported/confirmed
% 2 = implied by hospitals
cases_type = 1;

%% inputs
overlay = 31;
dateFrom = dd(2020,9,1);
dateTo = dd(2020,12,31);
data = struct();
data.rho = resize(db_age.Old/db_age.Total,dateFrom:dateTo);
data.X_obs = resize(db.NewCases,dateFrom:dateTo);
data.AC = resize(db.ActiveCases,dateFrom:dateTo);
data.TC = resize(db.TotalCases,dateFrom:dateTo);

suff = {'_default','_implied'};

info = load_fanchart_tseries('src_dir', '../R/Rt_estimation/results',...
        'src_filenames', {'output_R_reported.csv','output_R_implied.csv'},...
        'tar_dir','results','tar_filenames', {'Rt_reported.csv','Rt_implied.csv'});
data.Rt = resize(info{cases_type}.mean_ts,dateFrom:dateTo);

%% run calculations

% make guess of initial values/train model
% p = make_init_guess(s,data,dateFrom,dateTo);
p = train_model(s,data,dateFrom,dateTo);

save(strcat('results/forecast_init',suff{cases_type},'.mat'),'p');

% ********** 1./ train model
validFrom = dateFrom;%+30;
validTo = dateTo;
q = get_fcast_init_data(p,validFrom);

% make simple cross-validation
p1=make_init_forecast(s,q,validFrom,validTo);

% check results
check_targets('dbs',{p,p1},'varlist',{'X_o_obs','X_y_obs','X_obs'},'titles',{'X_o_obs','X_y_obs','X_obs'},...
    'dateFrom',validFrom,'dateTo',validTo);
