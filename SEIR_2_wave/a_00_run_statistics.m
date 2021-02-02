%% cleanup
initialize;

%% definitions
filename_in = 'data/calibration/hosp_data.xlsx';
filename_out = 'results/optimal_fit.mat';

%% data processing & distributions estimation
tb = readtable(filename_in);
save('data/calibration/hosp_data_raw.mat','tb');
s = run_clinical_inputs_statistics(tb); 

%% saving stuff
save(filename_out,'s');
