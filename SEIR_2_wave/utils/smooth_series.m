function res = smooth_series(Y,varargin)

if isa(Y,'tseries')
    convert = true;
    Y0 = 0*Y;
    Y = double(Y);
else
    convert = false;
end
if isempty(varargin)
    w = 7;
else
    w = varargin{1};
end
if length(varargin)<=1
    type = 3;
else
    type = varargin{2};
end
if length(varargin)<=2
    ends = 1;
else
    ends = varargin{3};
end
if length(varargin)<=3
    indiff = true;
else
    indiff = varargin{4};
end

if ~indiff
    YY = cumsum(Y);
    y = smooth_series_inner(YY,w,type,ends);
    y0 = [y(1); y(1:end-1)]; y1 = y;
    res = y1-y0;
else
    res = smooth_series_inner(Y,w,type,ends);
end

if convert
    res = Y0+res;
end
end