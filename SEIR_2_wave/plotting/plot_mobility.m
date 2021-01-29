function []=plot_mobility(mob,dateFrom,dateTo)

threshold = 100;

% raw
figure('Name','Mobility');
fn = fieldnames(mob);
for i=1:length(fn)
    yy = interp(resize(mob.(fn{i}),dateFrom:dateTo),dateFrom:dateTo);
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

% smooth
figure('Name','Mobility, smooth data');
for i=1:length(fn)
    yy = interp(resize(mob.(fn{i}),dateFrom:dateTo),dateFrom:dateTo);
    zz = smooth_series(mov_median(yy));
    if i==length(fn)
        plot(zz,'linewidth',2,'Color','k');hold on;
    else
        plot(zz,'linewidth',1);hold on;
    end
end
plot(bench,'color',[0.15 0.15 0.15],'linestyle','--');
grid on;
ylabel('%');
title('Mobility, smooth data');
legend(fn);

end
