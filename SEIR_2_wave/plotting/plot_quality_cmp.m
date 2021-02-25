function [] = plot_quality_cmp(db,db_rep,dateFrom,dateTo)

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
pp1 = plot(100*resize(db.rho_real_smooth,dateFrom:dateTo),'linewidth',2);
pp2 = plot(100*resize(db.rho,dateFrom:dateTo),'linewidth',2);
grid on;
legend([pp1 pp2],{'Hospitals-implied','Reported (PCR tests)'});
ylabel('% of total cases');
title('Share of 65+ in new cases');
subplot(2,1,2);
pp1 = plot(100*resize(db.sympt_share_real,dateFrom:dateTo),'linewidth',2);
pp2 = plot(100*resize(db.sympt_share_ideal,dateFrom:dateTo),'linewidth',2);
grid on;
legend([pp1 pp2],{'Observed','Ideal'});
ylabel('% of total cases');
title('Share symptomatic new cases');


end