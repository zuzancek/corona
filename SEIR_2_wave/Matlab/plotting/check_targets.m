function []=check_targets(varargin)

ip = inputParser;
addParamValue(ip, 'dbs', {}, @iscell);%#ok<*NVREPL>
addParamValue(ip, 'varlist',{}, @iscell);
addParamValue(ip, 'titles',{}, @iscell);
addParamValue(ip, 'legend',{}, @iscell);
addParamValue(ip, 'dateFrom',dd(2020,9,1), @isnumeric);
addParamValue(ip, 'dateTo',dd(2020,12,31), @isnumeric);
parse(ip, varargin{:});
results = ip.Results;

dateFrom = results.dateFrom;
dateTo = results.dateTo;
tit = results.titles;
dbs = results.dbs;
varlist = results.varlist;
leg = results.legend;
n = length(dbs);
m = length(varlist);

for i=1:m
    var = varlist{i};
    figure;
    for j=1:n
    	plot(resize(dbs{j}.(var),dateFrom:dateTo),'linewidth',1); hold on;
    end
    grid on;
    legend(leg);
    title(tit{i});
end
    
end
