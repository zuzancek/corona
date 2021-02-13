%% cleanup
initialize;

run_hosp = false;
out_db = cell(2,1);

%% definitions
for i=1:2
    if i==1
        filename_in = 'data/calibration/hosp_data.xlsx';
        filename_out = 'results/optimal_fit_hosp.mat';
        filename_mat = 'data/calibration/hosp_data_raw.mat';
    else
        filename_in = 'data/calibration/icu_data.xlsx';
        filename_out = 'results/optimal_fit_icu.mat';
        filename_mat = 'data/calibration/icu_data_raw.mat';
    end

    %% data processing & distributions estimation
    tb = readtable(filename_in);
    save(filename_mat,'tb');
    s = run_clinical_inputs_statistics(tb); 

    %% saving stuff
    save(filename_out,'s');
    out_db{i} = s;
end

par = setparam();
s = run_clinical_statistics_mild(out_db,par);

