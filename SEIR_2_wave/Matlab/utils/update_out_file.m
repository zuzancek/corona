function [] = update_out_file(varargin)

ip = inputParser;
addParamValue(ip, 'outdirs', {'../R/cases_forecast/data','../R/Rt_estimation/data'}, @ischar);%#ok<*NVREPL>
addParamValue(ip, 'filename_tar', {'cases_mm.csv','cases_raw.csv','cases_smooth.csv','cases_implied_smooth.csv','cases_implied_raw.csv'}, @iscell);%#ok<*NVREPL>
addParamValue(ip, 'filename_src', 'cases.csv', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'data_src', {'results/inputs.mat','results/inputs.mat','results/inputs.mat','results/results_impl.mat','results/results_impl.mat'}, @iscell);%#ok<*NVREPL>
addParamValue(ip, 'var', {'cases_data.cases_pcr_mm','cases_data.cases_pcr','cases_data.cases_pcr_smooth','cases_implied_data.X_smooth_all','cases_implied_data.X_all'}, @iscell);%#ok<*NVREPL>
addParamValue(ip, 'reported',true, @islogical);
parse(ip, varargin{:});

results = ip.Results;
outdirs = results.outdirs;
filename_tar = results.filename_tar;
data_src = results.data_src;
filename_src = results.filename_src;
var = results.var;

for i=1:length(outdirs)
    filepath_src = fullfile(pwd,outdirs{i},filename_src);
    T = readtable(filepath_src,'Delimiter',';','ReadVariableNames',true);
    src_dates = T.Datum;
    for j=1:length(data_src)        
        src = load(data_src{j});
        v = split(var{j},'.'); x = src;
        for k=1:length(v)
            x = x.(v{k});
        end
        xdates = src.dates;
        t1_idx = find(datetime(xdates.t1,'ConvertFrom','datenum')==src_dates);
        t2_idx = find(datetime(xdates.t2,'ConvertFrom','datenum')==src_dates);
        T.Dennych_PCR_prirastkov(t1_idx:t2_idx) = x(xdates.t1:xdates.t2);
        T.Pocet_potvrdenych_PCR_testami = cumsum(T.Dennych_PCR_prirastkov);
        T(t2_idx+1:end,:) = [];
        filepath_tar = fullfile(pwd,outdirs{i},filename_tar{j});
        writetable(T,filepath_tar,'Delimiter',';');
    end
end