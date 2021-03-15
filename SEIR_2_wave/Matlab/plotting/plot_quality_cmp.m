function [] = plot_quality_cmp(db,db_rep,s,dateFrom,dateTo) %#ok<INUSL>

%%
figure('Name','Testing quality I.');
bar(100*resize(db.obs_ratio_adj_raw,dateFrom:dateTo),'FaceAlpha',0.85);hold on;
pp1 = plot(100*resize(db.obs_ratio_adj,dateFrom:dateTo),'b','linewidth',2);
pp2 = plot(100*resize(db.obs_ratio_ideal,dateFrom:dateTo),'r','linewidth',2);
grid on;
legend([pp1 pp2],{'Real','Ideal'});
ylabel('% of total cases');
title('Observable ratio: real vs. ideal');

%%
figure('Name','Testing quality II.');
subplot(2,1,1);
bar(100*resize(db.X_o_smooth./db.X_smooth,dateFrom:dateTo),'FaceAlpha',0.75);hold on;
bar(100*resize(db.X_o_rep_smooth./db.X_rep_smooth,dateFrom:dateTo),'FaceAlpha',0.75);hold on;
pp1 = plot(100*resize(smooth_series(db.X_o_smooth./db.X_smooth),dateFrom:dateTo),'b','linewidth',2); hold on;
pp2 = plot(100*resize(smooth_series(db.X_o_rep_smooth./db.X_rep_smooth),dateFrom:dateTo),'r','linewidth',2);
pp3 = plot(100*(s.dep_ratio_65+0*resize(smooth_series(db.X_o_rep_smooth./db.X_rep_smooth),dateFrom:dateTo)),'k','linewidth',2);
grid on;
legend([pp1 pp2 pp3],{'Hospitals-implied','Reported (PCR tests)','Share in population'});
ylabel('% of total cases');
title('Share of 65+ in new cases');
subplot(2,1,2);
pp1 = plot(100*resize(db.sympt_share_real,dateFrom:dateTo),'linewidth',2); hold on;
pp2 = plot(100*resize(db.sympt_share_ideal,dateFrom:dateTo),'linewidth',2);
grid on;
legend([pp1 pp2],{'Observed','Ideal'});
ylabel('% of total cases');
title('Share symptomatic new cases');

%% 
figure('Name','Testing quality III.');
subplot(2,1,1);
bar(db.X_y_raw);hold on;
bar(db.X_o_raw);hold on;
pp1=plot(db.X_y_smooth,'b','linewidth',2);
pp2=plot(db.X_o_smooth,'r','linewidth',2);
pp3=plot(db.X_y_rep_smooth,'c-.','linewidth',2);
pp4=plot(db.X_o_rep_smooth,'m-.','linewidth',2);
grid on;
legend([pp1 pp2 pp3 pp4],{'Young - implied','Old - implied','Young - reported','Old - reported'});
ylabel('new cases');
title('Age structure of new cases: reported vs implied');

subplot(2,1,2);
bar(db.X_y_raw-db.X_y_rep_raw);hold on;
bar(db.X_o_raw-db.X_o_rep_raw);hold on;
pp1=plot(db.X_y_smooth-db.X_y_rep_smooth,'b','linewidth',2);
pp2=plot(db.X_o_smooth-db.X_o_rep_smooth,'r','linewidth',2);
grid on;
legend([pp1 pp2],{'Young','Old'});
ylabel('loast cases');
title('Lost cases (implied-reported)');

end