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

figure('Name','Clinical statistics');
%
subplot(2,2,1)
if mm
    bar(resize(data.H,dateFrom:dateTo));hold on;
end
if smooth
    plot(resize(data.H_smooth,dateFrom:dateTo),'linewidth',2);hold on;
end
if raw
    plot(resize(data.H_raw,dateFrom:dateTo),'linestyle','-.','Color',[0.5 0.5 0.5]);hold on;
end
grid on;
title('Total hospitalizations');
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
end
if smooth
    plot(resize(data.D_smooth,dateFrom:dateTo),'linewidth',2);hold on;
end
if raw
    plot(resize(data.D_raw,dateFrom:dateTo),'linestyle','-.','Color',[0.5 0.5 0.5]);hold on;
end
grid on;
title('Cummulative deaths (on+with)');

end