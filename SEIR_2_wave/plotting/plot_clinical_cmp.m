function [] = plot_clinical_cmp(db,db_alt,db_rep,dateFrom,dateTo,varargin)

ip = inputParser;
addParamValue(ip, 'mm', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'reported',true, @islogical);
parse(ip, varargin{:});
results = ip.Results;
mm = results.mm;
reported = results.reported;

%%
figure('Name','Situation in Hospitals I.');
subplot(2,1,1);
if mm
    bar(resize(db.H,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    bar(resize(db_alt.H,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
end
pp1=plot(resize(smooth_series(db.H),dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
pp2=plot(resize(smooth_series(db_alt.H),dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
if reported
    pp3=plot(resize(db_rep.H_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
    legend([pp1,pp2,pp3],{'reconstructed from implied cases','implied by reported daily new cases','observed'});
else
    legend([pp1,pp2],{'reconstructed from implied cases','implied by reported daily new cases'});
end
grid on;
title('Hospitalisations (total)');  

subplot(2,1,2);
if mm
    bar(resize(db.D,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    bar(resize(db_alt.D,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
end
pp1=plot(resize(smooth_series(db.D),dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
pp2=plot(resize(smooth_series(db_alt.D),dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
if reported
    pp3=plot(resize(db_rep.D_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
    legend([pp1,pp2,pp3],{'reconstructed from implied cases','implied by reported daily new cases','observed'});
else
    legend([pp1,pp2],{'reconstructed from implied cases','implied by reported daily new cases'});
end
grid on;
title('Deaths, cummulative (on+with)');  

%%
figure('Name','Situation in Hospitals II.');
subplot(2,1,1);
if mm
    bar(resize(db.S,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    bar(resize(db_alt.S,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
end
pp1=plot(resize(smooth_series(db.S),dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
pp2=plot(resize(smooth_series(db_alt.S),dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
if reported
    pp3=plot(resize(db_rep.S_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
    legend([pp1,pp2,pp3],{'reconstructed from implied cases','implied by reported daily new cases','observed'});
else
    legend([pp1,pp2],{'reconstructed from implied cases','implied by reported daily new cases'});
end
grid on;
title('Serious conditions (ICU,ECMO,V)');  

subplot(2,1,2);
if mm
    bar(resize(db.M,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
    bar(resize(db_alt.M,dateFrom:dateTo),'FaceAlpha',0.5,'EdgeAlpha',0.5);hold on;
end
pp1=plot(resize(smooth_series(db.M),dateFrom:dateTo),'linewidth',2,'Color','b');hold on;
pp2=plot(resize(smooth_series(db_alt.M),dateFrom:dateTo),'linewidth',2,'Color','r');hold on;
if reported
    pp3=plot(resize(db_rep.H_smooth-db_rep.S_smooth,dateFrom:dateTo),'linewidth',1,'Color','k','linestyle','--');hold on;
    legend([pp1,pp2,pp3],{'reconstructed from implied cases','implied by reported daily new cases','observed'});
else
    legend([pp1,pp2],{'reconstructed from implied cases','implied by reported daily new cases'});
end
grid on;
title('Moderate conditions (no ICU,ECMO,V needed)'); 

%%
figure('Name','Situation in Hospitals III.');
subplot(2,1,1)
pp1=bar(resize(db_rep.A,dateFrom:dateTo),'FaceAlpha',0.5);hold on;
pp2=plot(resize(db.IH,dateFrom:dateTo),'linewidth',2);
pp3=plot(resize(db_alt.IH,dateFrom:dateTo),'linewidth',2,'Color',[1 1 1]*0.5);
grid on;
legend([pp1 pp2 pp3],{'Reported','reconstructed from implied cases','implied by reported daily new cases'});
title('Hospital Admissions');

subplot(2,1,2)
pp1=bar(resize(db_rep.R,dateFrom:dateTo),'FaceAlpha',0.5);hold on;
pp2=plot(resize(db.HR,dateFrom:dateTo),'linewidth',2);
pp3=plot(resize(db_alt.HR,dateFrom:dateTo),'linewidth',2,'Color',[1 1 1]*0.5);
grid on;
legend([pp1 pp2 pp3],{'Reported','reconstructed from implied cases','implied by reported daily new cases'});
title('Hospital Discharges');

end