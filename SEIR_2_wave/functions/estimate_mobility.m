function estimate_mobility(mob,fcast_days,t0,t1,tw2)

threshold = 100;
s = setparam();

figure;
fn = fieldnames(mob);
for i=1:length(fn)
    yy = interp(resize(mob.(fn{i}),t0:t1),t0:t1);
    if i==length(fn)
        plot(yy,'linewidth',2,'Color','k');hold on;
    else
        plot(yy,'linewidth',1);hold on;
    end
end
xl = xlim;
bench = tseries(xl(1):xl(2),0)+threshold;
plot(bench,'color',[0.15 0.15 0.15],'linestyle','--');
grid on;
ylabel('%');
title('Mobility, raw data');
legend(fn);

figure;
for i=1:length(fn)
    yy = interp(resize(mob.(fn{i}),t0:t1),t0:t1);
    zz = smooth_series(double(yy),s.smooth_width,s.smooth_type,s.smooth_ends);
    yy = 0*yy+zz;
    if i==length(fn)
        plot(yy,'linewidth',2,'Color','k');hold on;
    else
        plot(yy,'linewidth',1);hold on;
    end
end
plot(bench,'color',[0.15 0.15 0.15],'linestyle','--');
grid on;
ylabel('%');
title('Mobility, smooth data');
legend(fn);

% w1 minimum
w1 = resize(yy,t0:tw2);
[min1,dw1] = min(w1);
thd1 = min(find(resize(w1,dw1:tw2)>threshold));
dt1 = thd1-dw1;
w10 = w1(dw1:thd1-1);w11 = w1(dw1+1:thd1);
gr1 = w11./w10-1;
gr1_m = 100*mean(gr1);
delta1 = (threshold-min1)/dt1;

% w2 minimum
w2 = resize(yy,tw2:t1);
tt1 = enddate(w2);
[min2,dw2] = min(w2);
w20 = w2(dw2:tt1-1);w21 = w2(dw2+1:tt1);
gr2 = w21./w20-1;
gr2_m = 100*mean(gr2);
delta2 = log(w2(tt1)-min2)/log(tt1-dw2);

figure;
plot(w11,'linewidth',1);hold on;
plot(w21,'linewidth',1);hold on;
legend({'1st wave','2nd wave'});
grid on;
title('Mobility (aggregate level): from minimum to 100%');

z0 = min2*(gr2_m/100+1)^(tt1-dw2+1+fcast_days);
y0 = tseries(tw2:tt1+fcast_days,0);
y0(tw2:tt1) = w2; y0(end) = z0; y0(tt1+1:end-1) = NaN; y0 = interp(y0,tw2:tt1+fcast_days);
z1 = min2*(0.5*(gr2(end)+gr2(end-1))+1)^(tt1-dw2+1+fcast_days);
y1 = tseries(tw2:tt1+fcast_days,0);
y1(tw2:tt1) = w2; y1(end) = z1; y1(tt1+1:end-1) = NaN; y1 = interp(y1,tw2:tt1+fcast_days,'method','pchip');
y2 = tseries(tw2:tt1+fcast_days,0);
y2(dw2:end) = (0:tt1-dw2+fcast_days).^delta2;  
y2(tw2:tt1) = w2; 
y2(tt1+1:end) = w2(dw2)+y2(tt1+1:end).*(1-0*[0:fcast_days-1]'/fcast_days);
y3 = tseries(tw2:tt1+fcast_days,0);
y3(dw2:end) = (0:tt1-dw2+fcast_days).^delta2;  
y3(tw2:tt1) = w2; 
y3(tt1+1:end) = w2(dw2)+y3(tt1+1:end).*(1-0.15*[0:fcast_days-1]'/fcast_days);

figure;
plot(y0);hold on; 
plot(y1);hold on; 
plot(y2);hold on; 
plot(y3);hold on; 

scatter(tt1,y0(tt1),'k');
grid on;

end