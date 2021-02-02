%% cleanup
initialize;

%% definitions
filename_in = 'data/calibration/hosp_data.xlsx';
filename_out = 'results/optimal_fit.mat';

%% data processing & distributions estimation
s = run_clinical_inputs_statistics(filename_in); 

%% saving stuff
save(filename_out,'s');
