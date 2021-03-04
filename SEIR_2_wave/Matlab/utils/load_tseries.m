function []=load_tseries(varargin)

ip = inputParser;
addParamValue(ip, 'src_dir', {'../R/Rt_estimation/results/'}, @ischar);%#ok<*NVREPL>
addParamValue(ip, 'src_filenames', {'cases_mm.csv','cases_implied.csv'}, @iscell);%#ok<*NVREPL>
addParamValue(ip, 'tar_dir', 'results/', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'tar_filenames', {'Rt_reported.csv','Rt_implied.csv'}, @ischar);%#ok<*NVREPL>
parse(ip, varargin{:});