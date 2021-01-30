function []=plot_epidemiology_reporting(cases_data,test_data,dateFrom,dateTo,varargin)

ip = inputParser;
addParamValue(ip, 'raw', false, @islogical); %#ok<*NVREPL>
addParamValue(ip, 'mm', true, @islogical);
addParamValue(ip, 'smooth',true, @islogical);
parse(ip, varargin{:});
results = ip.Results;
raw = results.raw;
mm = results.mm;
smooth = results.smooth;

cases_data.cases_ag = dI_inflow_ag;   cases_data.cases_ag_mm = mov_median(dI_inflow_ag);   cases_data.cases_ag_smooth = smooth_series(cases_data.cases_ag_mm);
cases_data.cases_pcr = dI_inflow_pcr; cases_data.cases_pcr_mm = mov_median(dI_inflow_pcr); cases_data.cases_pcr_smooth = smooth_series(cases_data.cases_pcr_mm);
cases_data.cases_total = dI_inflow;   cases_data.cases_total_mm = mov_median(dI_inflow);   cases_data.cases_total_smooth = smooth_series(cases_data.cases_total_mm);
cases_data.ptr_pcr = pos_test_ratio;  cases_data.ptr_mm = mov_median(pos_test_ratio);      cases_data.ptr_smooth = smooth_series(cases_data.ptr_mm);
cases_data.ptr_ag = pos_test_ratio_ag;cases_data.ptr_ag_mm = mov_median(pos_test_ratio_ag);cases_data.ptr_ag_smooth = smooth_series(cases_data.ptr_ag_mm);

%
figure('Name','New confirmed cases, PCR, AG')
leg = {};
if mm
    bar(resize(cases_data.cases_total,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    bar(resize(cases_data.cases_pcr,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    leg = {leg{:},{'PCR only (7d median)','PCR+AG (7d median)'}}; %#ok<*CCAT>
end
if smooth
    plot(resize(cases_data.cases_total_smooth,dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
    plot(resize(cases_data.cases_pcr_smooth,dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
    leg = {leg{:},{'PCR only (smooth)','PCR+AG (smooth)'}};
end
if raw
    plot(resize(cases_data.cases_total_raw,dateFrom:dateTo),'linewidth',1,'Color',[0 0.25 0.75],'linestyle','-.');hold on;
    plot(resize(cases_data.cases_pcr_raw,dateFrom:dateTo),'linewidth',1,'Color',[0.75 0.25 0],'linestyle','-.');hold on;
    leg = {leg{:},{'PCR only (raw)','PCR+AG (raw)'}};
end
title('New infections (confirmed only)');
legend(leg);
grid on;

figure('Name','Testing characteristics I.');
subplot(2,1,1)
plot(resize(test_data.ptr_pcr,dateFrom:dateTo),'linewidth',1); hold on;
plot(resize(test_data.ptr_ag,dateFrom:dateTo),'linewidth',1); hold on;
pp1=plot(resize(test_data.ptr_pcr_smooth,dateFrom:dateTo),'Color','b','linewidth',2);
pp2=plot(resize(test_data.ptr_ag_smooth,dateFrom:dateTo),'Color','r','linewidth',2);
title('Positive tests ratio');
legend([pp1 pp2],{'PCR','AG'});
grid on;
subplot(2,1,2);
plot(resize(test_data.tests_pcr,dateFrom:dateTo),'linewidth',1);hold on;
plot(resize(test_data.tests_ag,dateFrom:dateTo),'linewidth',1);hold on;
pp1=plot(resize(test_data.tests_pcr_smooth,dateFrom:dateTo),'Color','b','linewidth',2);
pp2=plot(resize(test_data.tests_ag_smooth,dateFrom:dateTo),'Color','r','linewidth',2);
title('Tests (reported)');
legend([pp1 pp2],{'PCR','AG'});
grid on;

%
figure('Name','Testing characteristics II.')
subplot(2,1,1);
plot(resize(100*cases_data.asymp_ratio,dateFrom:dateTo),'linewidth',1);hold on;
plot(resize(100*cases_data.asymp_ratio_smooth,dateFrom:dateTo),'linewidth',2);hold on;
title('Share of asymptomatic new cases (confirmed PCR only)');
legend({'raw','smooth'}); 
ylabel('%');
grid on;
subplot(2,1,2);
plot(resize(100*cases_data.old_ratio,dateFrom:dateTo),'linewidth',1);hold on;
plot(resize(100*cases_data.old_ratio_smooth,dateFrom:dateTo),'linewidth',2);hold on;
title('Share of 65+ positive cases (confirmed PCR only)');
legend({'raw','smooth'}); 
ylabel('%');
grid on;

end