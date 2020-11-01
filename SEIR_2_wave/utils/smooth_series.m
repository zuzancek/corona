function res = smooth_series(Y,w,type,ends,varargin)

try
    indiff = (varargin{1} == false);
catch err %#ok<NASGU>
    indiff = true;
end

if ~indiff
    YY = cumsum(Y);
    y = smooth_series_inner(YY,w,type,ends);
    y0 = [y(1); y(1:end-1)]; y1 = y;
    res = y1-y0;
else
    res = smooth_series_inner(Y,w,type,ends);
end

end