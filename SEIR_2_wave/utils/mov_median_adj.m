function [y] = mov_median_adj(x,varargin)

try
    d = varargin{1};
    try
        dateFrom = varargin{2};
        try
            dateTo = varargin{3};
        catch err2 %#ok<*NASGU>
            if isa(x,'tseries')
                dateTo = enddate(x);
            else
                dateTo = length(x);
            end
        end
    catch err1
        if isa(x,'tseries')
            dateFrom = startdate(x);
            dateTo = enddate(x);
        else
            dateFrom = 1;
            dateTo = length(x);
        end
    end
catch err0
    d = [7 0];
    if isa(x,'tseries')
        dateFrom = startdate(x);
        dateTo = enddate(x);
    else
        dateFrom = 1;
        dateTo = length(x);
    end
end

y0 = movmedian(x(dateFrom:dateTo,:),d,'omitnan');
y1 = 0.5*(y0(floor(d/2):end-1,:)+y0(ceil(d/2):end,:));
y = smooth_series(x(dateFrom:dateTo,:),d(1));
y(1:length(y1)-1,:) = y1(2:end,:);

if isa(x,'tseries')
    y = y+0*x;
end

end

