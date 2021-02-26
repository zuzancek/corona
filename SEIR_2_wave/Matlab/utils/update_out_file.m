function [] = update_out_file(varargin)

ip = inputParser;
addParamValue(ip, 'outdir', '../R/cases_forecast/data', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'filename_tar', 'cases_implied.csv', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'filename_src', 'cases.csv', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'data_filename_src', 'results/results_impl.mat', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'reported',true, @islogical);
parse(ip, varargin{:});


results = ip.Results;
outdir = results.outdir;
filename_tar = results.filename_tar;
filename_src = results.filename_src;
data_filename_src = results.data_filename_src;

filepath_src = fullfile(pwd,outdir,filename_src);
T = readtable(filepath_src,'Delimiter',';','ReadVariableNames',true);

src_dates = T.Datum;
data_src = load(data_filename_src);
x = data_src.cases_implied_data.X_smooth_all;
xdates = data_src.dates;
t1_idx = find(datetime(xdates.t1,'ConvertFrom','datenum')==src_dates);
t2_idx = find(datetime(xdates.t2,'ConvertFrom','datenum')==src_dates);
T.Dennych_PCR_prirastkov(t1_idx:t2_idx) = x(xdates.t1:xdates.t2);
T.Pocet_potvrdenych_PCR_testami = cumsum(T.Dennych_PCR_prirastkov);
T(t2_idx:end,:) = [];

filepath_tar = fullfile(pwd,outdir,filename_tar);
writetable(T,filepath_tar,'Delimiter',';');

end