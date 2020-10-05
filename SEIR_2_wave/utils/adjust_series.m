function [x_smooth,x_new,y_new] = adjust_series(x_in,s_width,s_type,s_end,scale,threshold,y)

y_new = max(y,threshold);
x_new = ceil(x_in.*(1+scale.*(y>threshold)));
x_smooth = smooth_series(x_new,s_width,s_type,s_end);

end