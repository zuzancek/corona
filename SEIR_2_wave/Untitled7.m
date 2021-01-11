close all;
x = 0:20;
xx = 0:0.01:x(end);
d_std = 0.62; d_mean = 6.5;
d_shape = d_mean*d_std^2; d_scale = 1/d_std^2;
y = pdf('Gamma',x,d_shape,d_scale);
yy = pdf('Gamma',xx,d_shape,d_scale);
figure;
plot(xx,yy);hold on;
plot(x,y,'o')
s = sum(y);
y = y/s;
s = '';
for i=1:length(x)
    s = strcat(s,'%2.2f\t');
end
s = strcat(s,'\n');
fprintf(s,x);
fprintf(s,100*y);