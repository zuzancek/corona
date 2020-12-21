function [z,z_smooth,z_ext,z_ext_smooth] = process_xls(filename,dateFrom,dateTo,extendFrom,extendTo,s,initial,final)

tbl = readtable(filename);
x = tbl.x;
y = tbl.y;
xx = 1:dateTo-dateFrom+1;
x1 = ceil(max(x));
d = 1;
if x1==x(end)
    x(end) = [];
    d = 0;
end
x_new = [x;[x1:xx(end)]'];
y_new = [y;[y(end)+zeros(length(xx)-x1+d,1)]];
yy = interp1(x_new,y_new,xx,'pchip');
z = tseries(dateFrom:dateTo,yy);
z_ext = tseries(extendFrom:extendTo,NaN);
z_ext(dateFrom:dateTo) = z;
if isempty(initial)
    initial = z(dateFrom);
end
if isempty(final)
    final = z(dateTo);    
end

z_ext(extendFrom) = initial; z_ext(extendTo) = final;
z_ext = interp(z_ext,extendFrom:extendTo);
z_ext_smooth = smooth_series(z_ext,s.smooth_width,s.smooth_type,s.smooth_ends);
z_smooth = smooth_series(z,s.smooth_width,s.smooth_type,s.smooth_ends);

z = z/100;
z_ext = z_ext/100;
z_smooth = z_smooth/100;
z_ext_smooth = z_ext_smooth/100;

end