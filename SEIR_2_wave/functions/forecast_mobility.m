function [mob] = forecast_mobility(mob,fcast_days,t0,t1,tw2)

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
%
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

figure;
yy = interp(resize(mob.SK,t0:t1),t0:t1);
plot(yy,'linewidth',1);hold on;
for k=4:7
    zz = smooth_series(double(yy),k,s.smooth_type,s.smooth_ends);
    plot(0*yy+zz,'linewidth',1);hold on;
end
udata = 0*yy+zz;
plot(bench,'color',[0.15 0.15 0.15],'linestyle','--');
grid on;
legend({'raw','smooth (k=4)','smooth (k=5)','smooth (k=6)','smooth (k=7)'});
ylabel('%');
title('Mobility, aggregate level');

% w1 minimum
w1 = resize(yy,t0:tw2);
[min1,dw1] = min(w1);
thd1 = min(find(resize(w1,dw1:tw2)>threshold)); %#ok<*MXFND>
dt1 = thd1-dw1;
w10 = w1(dw1:thd1-1);w11 = w1(dw1+1:thd1);
gr1 = w11./w10-1;
gr1_m = 100*mean(gr1); %#ok<*NASGU>
delta1 = (threshold-min1)/dt1;

% w2 minimum
w2 = resize(yy,tw2:t1);
tt1 = enddate(w2);
[min2,dw2] = min(w2);
w20 = w2(dw2:tt1-1);w21 = w2(dw2+1:tt1);
gr2 = w21./w20-1;
gr2_m = 100*mean(gr2);
delta2 = log(w2(tt1)-min2)/log(tt1-dw2);

z0 = min2*(gr2_m/100+1)^(tt1-dw2+1+fcast_days);
y0 = tseries(tw2:tt1+fcast_days,0);
y0(tw2:tt1) = w2; y0(end) = z0; y0(tt1+1:end-1) = NaN; y0 = interp(y0,tw2:tt1+fcast_days);
zz = smooth_series(double(y0),s.smooth_width,s.smooth_type,s.smooth_ends);
y0 = 0*y0+zz;
z1 = min2*(0.5*(gr2(end)+gr2(end-1))+1)^(tt1-dw2+1+fcast_days);
y1 = tseries(tw2:tt1+fcast_days,0);
y1(tw2:tt1) = w2; y1(end) = z1; y1(tt1+1:end-1) = NaN; y1 = interp(y1,tw2:tt1+fcast_days,'method','pchip');
zz = smooth_series(double(y1),s.smooth_width,s.smooth_type,s.smooth_ends);
y1 = 0*y1+zz;
y2 = tseries(tw2:tt1+fcast_days,0);
y2(dw2:end) = (0:tt1-dw2+fcast_days).^delta2;  
y2(tw2:tt1) = w2; 
y2(tt1+1:end) = w2(dw2)+y2(tt1+1:end).*(1-0*[0:fcast_days-1]'/fcast_days); %#ok<*NBRAK>
zz = smooth_series(double(y2),s.smooth_width,s.smooth_type,s.smooth_ends);
y2 = 0*y2+zz;
y3 = tseries(tw2:tt1+fcast_days,0);
y3(dw2:end) = (0:tt1-dw2+fcast_days).^delta2;  
y3(tw2:tt1) = w2; 
y3(tt1+1:end) = w2(dw2)+y3(tt1+1:end).*(1-0.3*[0:fcast_days-1]'/fcast_days);
zz = smooth_series(double(y3),s.smooth_width,s.smooth_type,s.smooth_ends);
y3 = 0*y3+zz;
w2s = resize(y3,tw2:tt1);
y4 = tseries(tw2:tt1+fcast_days,0);
y4(tt1:end) = w2(tt1);  
y4(tw2:tt1) = w2; 
zz = smooth_series(double(y4),s.smooth_width,s.smooth_type,s.smooth_ends);
y4 = 0*y4+zz;

figure;
plot(y3,'Color',[0 0.4470 0.7410],'linewidth',1);hold on; 
ylim([50 150]);
yl = ylim();
area([tw2 tt1],[yl(2) yl(2)],'FaceColor',[0.75 0.75 0.75],'EdgeColor','k','FaceAlpha',0.25,'EdgeAlpha',0);
pp1 = plot(y3,'Color',[0 0.4470 0.7410],'linewidth',1);hold on; 
% plot(y2);hold on; 
pp2 = plot(y4,'linewidth',1,'Color',[0.8500 0.3250 0.0980]);hold on;
yy = 0.5*(y4+y3);
pp3 = plot(yy,'linewidth',1,'Color',[0.9290 0.6940 0.1250]);hold on;
plot(w2s,'linewidth',1,'Color',[0.5 0.5 0.5]);hold on; 
scatter(tt1,y0(tt1),'k*');
xl = xlim();
bench = tseries(xl(1):xl(2),0)+threshold;
plot(bench,'color',[0.15 0.15 0.15],'linestyle','--');
grid on;
ylabel('%');
title('Mobility (aggregate): forecast');
legend([pp1 pp2 pp3],{'High','Low','Medium'});

mob = struct();
mob.low = y3;
mob.medium = yy;
mob.high = y0;
mob.dateFrom = tt1+1;
mob.days = fcast_days;
mob.low_norm = y3/y3(tt1);
mob.medium_norm = yy/yy(tt1);
mob.high_norm = y0/y0(tt1);
mob.orig = udata;

end