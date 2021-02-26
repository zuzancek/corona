function [] = plot_clinical_statistics(data,dateFrom,dateTo,varargin)

ip = inputParser;
addParamValue(ip, 'raw', false, @islogical); %#ok<*NVREPL>
addParamValue(ip, 'mm', true, @islogical);
addParamValue(ip, 'smooth',true, @islogical);
parse(ip, varargin{:});
results = ip.Results;
raw = results.raw;
mm = results.mm;
smooth = results.smooth;

figure('Name','Clinical statistics I.');
%
subplot(2,2,1)
if mm
    bar(resize(data.H,dateFrom:dateTo));hold on;
    bar(resize(data.H_c,dateFrom:dateTo));hold on;
    bar(resize(data.H_s,dateFrom:dateTo),'FaceColor',[0.5 0.5 0.5]);hold on;
end
if smooth
    plot(resize(data.H_smooth,dateFrom:dateTo),'b','linewidth',2);hold on;
    plot(resize(data.H_c_smooth,dateFrom:dateTo),'r','linewidth',2);hold on;
    plot(resize(data.H_s_smooth,dateFrom:dateTo),'Color',0.75*[1 1 1],'linewidth',2);hold on;
end
grid on;
legend('total','confirmed','suspected');
title('Hospitalizations');
%
subplot(2,2,2)
if mm
    bar(resize(data.C,dateFrom:dateTo));hold on;
end
if smooth
    plot(resize(data.C_smooth,dateFrom:dateTo),'linewidth',2);hold on;
end
if raw
    plot(resize(data.C_raw,dateFrom:dateTo),'linestyle','-.','Color',[0.5 0.5 0.5]);hold on;
end
grid on;
title('ICU hospitalizations');
%
subplot(2,2,3)
if mm
    bar(resize(data.V,dateFrom:dateTo));hold on;
end
if smooth
    plot(resize(data.V_smooth,dateFrom:dateTo),'linewidth',2);hold on;
end
if raw
    plot(resize(data.V_raw,dateFrom:dateTo),'linestyle','-.','Color',[0.5 0.5 0.5]);hold on;
end
grid on;
title('Ventilations');
%
subplot(2,2,4)
if mm
    bar(resize(data.D,dateFrom:dateTo));hold on;
    bar(resize(data.D_on,dateFrom:dateTo));hold on;
    bar(resize(data.D_with,dateFrom:dateTo),'FaceColor',[0.5 0.5 0.5]);hold on;
end
if smooth
    plot(resize(data.D_smooth,dateFrom:dateTo),'b','linewidth',2);hold on;
    plot(resize(data.D_on_smooth,dateFrom:dateTo),'r','linewidth',2);hold on;
    plot(resize(data.D_with_smooth,dateFrom:dateTo),'Color',0.75*[1 1 1],'linewidth',2);hold on;
end
grid on;
legend('total','on','with');
title('Cummulative deaths');


figure('Name','Clinical statistics II.');
%
subplot(2,1,1)
if mm
    bar(100*resize(data.death_old_ratio,dateFrom:dateTo));hold on;
end
if smooth
    plot(100*resize(data.death_old_ratio_smooth,dateFrom:dateTo),'b','linewidth',2);hold on;
end
grid on;
ylabel('%');
title('Share of 65+ in deaths');

subplot(2,1,2)
if mm
    bar(100*resize(1-data.hosp_y_ratio,dateFrom:dateTo));hold on;
end
if smooth
    plot(100*smooth_series(resize(1-data.hosp_y_ratio,dateFrom:dateTo)),'b','linewidth',2);hold on;
end
ylabel('%');
grid on;
title('Share of 65+ in hospitals');

end