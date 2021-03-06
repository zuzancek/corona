%% cleanup
initialize;

run_hosp = false;
out_db = cell(2,1);
filename_out = 'results/optimal_fit.mat';

%% definitions

for i=1:2
    if i==1
        filename_in = 'data/calibration/hosp_data.xlsx';
        filename_mat = 'data/calibration/hosp_data_raw.mat';
    else
        filename_in = 'data/calibration/icu_data.xlsx';
        filename_mat = 'data/calibration/icu_data_raw.mat';
    end

    %% data processing & distributions estimation
    tb = readtable(filename_in);
    save(filename_mat,'tb');
    s = run_clinical_inputs_statistics(tb); 

    %% saving stuff
    out_db{i} = s;
end

%  calculate remaining statistics
db = process_clinical_inputs_statistics(out_db);
stat_total = db{1};  
stat_severe = out_db{2};
stat_mild = db{2}; stat_mild_severe = db{3}; 
raw = out_db;
save(filename_out,'stat_total','stat_severe','stat_mild','stat_mild_severe','raw');
