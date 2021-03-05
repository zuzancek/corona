function [info] = load_fanchart_tseries(varargin)

ip = inputParser;
addParamValue(ip, 'src_dir', '../R/Rt_estimation/results', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'src_filenames', {'output_R_reported.csv','output_R_implied.csv'}, @iscell);%#ok<*NVREPL>
addParamValue(ip, 'tar_dir', 'results', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'tar_filenames', {'Rt_reported.csv','Rt_implied.csv'}, @iscell);%#ok<*NVREPL>
parse(ip, varargin{:});

results = ip.Results;
src_dir = results.src_dir;
src_filenames = results.src_filenames;
tar_dir = results.tar_dir;
tar_filenames = results.tar_filenames;
info = cell(length(src_filenames),1);

for i=1:length(src_filenames)
    filepath_src = fullfile(pwd,src_dir,src_filenames{i});
    filepath_tar = fullfile(pwd,tar_dir,tar_filenames{i});
    T = readtable(filepath_src,'Delimiter',',','ReadVariableNames',true);
    NT = size(T,1);
    namelist = T.Properties.VariableNames;
    idxcol = find(startsWith(namelist,'Q_'));   
    NQ = length(idxcol);
    info{i}.dateFrom = datenum(T.Date(1));
    info{i}.dateTo = datenum(T.Date(end));
    Q_mat = zeros(NT,NQ);
    Q_mat_ts = cell(NQ,1);
    for j=1:length(idxcol)
        vals=cellfun(@(x) strrep(x,'NA',''),T{:,idxcol(j)},'UniformOutput',false);
        vals=cellfun(@(x) str2double(x),vals);
        Q_mat(:,j) = vals;        
        Q_mat_ts{j} = tseries(info{i}.dateFrom:info{i}.dateTo,vals);
        db.(namelist{idxcol(j)}) = Q_mat_ts{j};
    end
    info{i}.Q_mat = Q_mat;
    info{i}.Q_mat_ts = Q_mat_ts;
    info{i}.idx0 = find(~isnan(Q_mat(:,1)),1);
    info{i}.X = T.Cases;
    info{i}.X_ts = tseries(info{i}.dateFrom:info{i}.dateTo,T.Cases);
    idxcol = [find(startsWith(namelist,'Median_')),find(startsWith(namelist,'Mean')),...
        find(startsWith(namelist,'Std_')),find(startsWith(namelist,'Cases_'))];   
    varlist = {'median','mean','std','X'};
    info{i}.dateFrom = datenum(T.Date(1));
    info{i}.dateTo = datenum(T.Date(end));
    for j=1:length(idxcol)
        vals=cellfun(@(x) strrep(x,'NA',''),T{:,idxcol(j)},'UniformOutput',false);
        vals=cellfun(@(x) str2double(x),vals);
        vals_ts = tseries(info{i}.dateFrom:info{i}.dateTo,vals);
        info{i}.(varlist{j}) = vals;
        info{i}.(strcat(varlist{j},'_ts')) = vals_ts;
        db.(varlist{j}) = vals_ts;
    end
    dbsave(db,filepath_tar);
end