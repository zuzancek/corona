function [] = plot_epidemiology_cmp(db,db_alt,db_rep,dateFrom,dateTo,varargin)

ip = inputParser;
addParamValue(ip, 'implied', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'reported',true, @islogical);
parse(ip, varargin{:});
results = ip.Results;
implied = results.implied;
reported = results.reported;

%%
figure('Name','Daily cases: reported vs.real');
pp1=bar(resize(db.X_rep_raw,dateFrom:dateTo),'FaceAlpha',0.85);hold on;
p=bar(resize(db.X_raw,dateFrom:dateTo-db.fcast_per-1),'FaceAlpha',0.7);
bar(resize(db.X_all,dateTo-db.fcast_per:dateTo),'FaceAlpha',0.4,'FaceColor',p.FaceColor);
pp2=bar(mov_median(resize(db_rep.H,dateFrom:dateTo)),'FaceAlpha',0.5,'FaceColor',[0.5 0.5 0.5]);
pp3=bar(mov_median(resize(db_rep.S,dateFrom:dateTo)),'FaceAlpha',0.5,'FaceColor','k');
plot(resize(db.X_rep_smooth,dateFrom:dateTo),'b','linewidth',2);
plot(resize(db.X_smooth,dateFrom:dateTo-db.fcast_per-1),'r','linewidth',2);
plot(resize(smooth_series(db.X_all),dateTo-db.fcast_per:dateTo),'r-.','linewidth',2);
plot(smooth_series(mov_median(resize(db_rep.H,dateFrom:dateTo))),'--','Color',[0.25 0.25 0.25],'linewidth',1);
plot(smooth_series(mov_median(resize(db_rep.S,dateFrom:dateTo))),'k-.','linewidth',1);
grid on;
legend([pp1 p pp2 pp3],{'New cases: officially reported','New cases: implied by hospitals', 'Hospitalizations', 'Intensive Care'});
title('New cases: reported vs. real');

%%
dt = dateTo-dateFrom;
if implied
    figure('Name','Hospitalization/Deaths: backcheck');
    pp3=bar(tseries(dateFrom:dateTo,db.H(end-dt:end))); hold on;
    pp5=bar(tseries(dateFrom:dateTo,db.D(end-dt:end)),'FaceAlpha',0.5);hold on;
    pp7=bar(tseries(dateFrom:dateTo,db.S(end-dt:end)),'FaceColor',0.25*[1 1 1]); hold on;
    pp1=plot(resize(smooth_series(db.X),dateFrom:dateTo),'linewidth',2,'Color',[0.75 0.75 0]);hold on; 
    pp2=plot((resize(db_rep.H,dateFrom:dateTo)),'b','linewidth',2);hold on;
    pp6=plot((resize(db_rep.S,dateFrom:dateTo)),'k','linewidth',2);hold on;
    pp4=plot((resize(db_rep.D,dateFrom:dateTo)),'r','linewidth',2);hold on;
    grid on;
    legend([pp1 pp2 pp3 pp4 pp5 pp6 pp7],{'New cases','Hospitalizations - reported', ...
        'Hospitalizations - implied', 'Deaths - reported, cummulative','Deaths - implied, cummulative', ...
        'Serious cases - reported', 'Serious cases - implied'});
end

if reported
    figure('Name','Hospitalizations/Deaths implied by reported data');
    pp3=bar(tseries(dateFrom:dateTo,db_alt.H(end-dt:end))); hold on;
    pp5=bar(tseries(dateFrom:dateTo,db_alt.D(end-dt:end)),'FaceAlpha',0.5);hold on;
    pp7=bar(tseries(dateFrom:dateTo,db_alt.S(end-dt:end)),'FaceColor',0.25*[1 1 1]); hold on;
    pp1=plot(smooth_series(resize(db_alt.X,dateFrom:dateTo)),'linewidth',2,'Color',[0.75 0.75 0]);hold on; 
    pp2=plot((resize(db_rep.H,dateFrom:dateTo)),'b','linewidth',2);hold on;
    pp6=plot((resize(db_rep.S,dateFrom:dateTo)),'k','linewidth',2);hold on;
    pp4=plot((resize(db_rep.D,dateFrom:dateTo)),'r','linewidth',2);hold on;
    grid on;
    legend([pp1 pp2 pp3 pp4 pp5 pp6 pp7],{'New cases','Hospitalizations - reported', ...
        'Hospitalizations - implied', 'Deaths - reported, cummulative','Deaths - implied, cummulative', ...
        'Serious cases - reported', 'Serious cases - implied'});
end

end