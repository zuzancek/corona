function [y] = mov_median(x,varargin)

try
    d = varargin{1};
    try
        dateFrom = varargin{2};
        try
            dateTo = varargin{3};
        catch err2
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

y = movmedian(x(dateFrom:dateTo,:),d,'omitnan');

if isa(x,'tseries')
    y = y+0*x;
end

end

