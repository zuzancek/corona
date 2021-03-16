
%% initialization
initialize;

%% load data
db = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
db_age = dbload('data/new_cases_age.csv','dateFormat','yyyy-mm-dd','freq','daily');