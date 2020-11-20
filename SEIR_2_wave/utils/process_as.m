function [z,z_smooth] = process_as(filename,dateFrom,dateTo,s,final)

tbl = readtable(filename);
x = tbl.x;
y = tbl.y;
xx = 1:dateTo-dateFrom+1;
yy = interp1(x,y,xx,'cubic');
z = tseries(dateFrom:dateTo,yy);
z_smooth = smooth_series(z,s.smooth_width,s.smooth_type,s.smooth_ends);

try 
    zz = tseries(dateFrom:final.date,NaN); 
    zz(dateFrom:dateTo) = z;
    zz(final.date) = final.value;
    zz = interp(zz);
    z = zz;
    z_smooth = smooth_series(z,s.smooth_width,s.smooth_type,s.smooth_ends);
catch err %#ok<NASGU>
end

z = z/100;
z_smooth = z_smooth/100;

end