function [y,y_smooth] = extend_series(x,t0,t1,v0,v1)

y = tseries(t0:t1,NaN);
tt0 = startdate(x);
tt1 = enddate(x);
y(tt0:tt1) = x;
if isempty(v0)
    v0 = x(tt0);
end
if isempty(v1)
    v1 = x(tt1);
end
y(t0) = v0;
y(t1) = v1;
y = interp(y,t0:t1);
y_smooth = smooth_series(y,7,5,1);

end