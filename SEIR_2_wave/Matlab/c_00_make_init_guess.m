
%% initialization
initialize;

%% load data
db = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_age = dbload('data/new_cases_age.csv','dateFormat','yyyy-mm-dd','freq','daily');
s = setparam();

%% inputs
dateFrom = dd(2020,8,1);
dateTo = dd(2020,12,31);
data = struct();
data.rho = resize(db_age.Old/db_age.Total,dateFrom:dateTo);
data.X_obs = resize(db.NewCases,dateFrom:dateTo);
data.AC = resize(db.ActiveCases,dateFrom:dateTo);
data.TC = resize(db.TotalCases,dateFrom:dateTo);

info = load_fanchart_tseries('src_dir', '../R/Rt_estimation/results',...
        'src_filenames', {'output_R_reported.csv','output_R_implied.csv'},...
        'tar_dir','results','tar_filenames', {'Rt_reported.csv','Rt_implied.csv'});
data.Rt = resize(info{1}.mean_ts,dateFrom:dateTo);

%% run calcaulations
p = make_init_guess(s,data,dateFrom,dateTo);