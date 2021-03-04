function []=load_fanchart_tseries(varargin)

ip = inputParser;
addParamValue(ip, 'src_dir', {'../R/Rt_estimation/results'}, @ischar);%#ok<*NVREPL>
addParamValue(ip, 'src_filenames', {'cases_mm.csv','cases_implied.csv'}, @iscell);%#ok<*NVREPL>
addParamValue(ip, 'tar_dir', 'results', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'tar_filenames', {'Rt_reported.csv','Rt_implied.csv'}, @ischar);%#ok<*NVREPL>
parse(ip, varargin{:});

results = ip.Results;
src_dir = results.src_dir;
src_filenames = results.src_filenames;
tar_dir = results.tar_dir;
tar_filenames = results.tar_filenames;
% var = {results.var_0,results.var_1};

% "","t_start","t_end","Mean(R)","Std(R)","Quantile.0.025(R)","Quantile.0.05(R)","Quantile.0.25(R)","Median(R)","Quantile.0.75(R)","Quantile.0.95(R)","Quantile.0.975(R)","Date","Cases"

for i=1:length(src_filenames)
    filepath_src = fullfile(src_dir,src_filenames{i});
    T = readtable(filepath_src,'Delimiter',',','ReadVariableNames',true);
    src_dates = T.Date;
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