initialize;
today = datetime(2021,05,24);

%% load
supply_astra=readtable('data/supply_plan_astrazeneca.csv');
supply_pfizer=readtable('data/supply_plan_pfizer.csv');
supply_jj=readtable('data/supply_plan_jj.csv');
supply_moderna=readtable('data/supply_plan_moderna.csv');
supply_total = sortrows(union(union(union(supply_astra,supply_pfizer,'stable'),supply_jj,'stable'),supply_moderna,'stable'),1);
vacc_agg = readtable('data/vacc_date.csv');
target = 3300000;

%% visualise
figure;
plot(supply_total.Date,cumsum(supply_total.Doses),'linewidth',1);hold on;
yyaxis right
bar(supply_astra.Date,(supply_astra.Doses),.5,'FaceColor',[.75 .25 0],'FaceAlpha',.5);
bar(supply_pfizer.Date,(supply_pfizer.Doses),'FaceColor',[.125 0.75 .125],'FaceAlpha',.5);
bar(supply_jj.Date,(supply_jj.Doses),.5,'FaceColor',[.25 0 0.5],'FaceAlpha',.5);
bar(supply_moderna.Date,(supply_moderna.Doses),.5,'FaceColor',[0 .5 .5],'FaceAlpha',.5);
grid on;title('Total supply');

figure;
subplot(2,2,1);
% firstdate = supply_astra.Date(1);
% cumsupply = cumsum(supply_astra.Doses);
% area(0:today-firstdate,cumsupply(end):cumsupply(end),'FaceColor',0.75*[1 1 1],'EdgeColor',0.65*[1 1 1],'FaceAlpha',0.5,'EdgeAlpha',0);hold on;
plot(supply_astra.Date,cumsum(supply_astra.Doses),'linewidth',1);hold on;
bar(supply_astra.Date,(supply_astra.Doses));
grid on;title('AstraZeneca: supply');
subplot(2,2,2);
plot(supply_pfizer.Date,cumsum(supply_pfizer.Doses),'linewidth',1);hold on;
bar(supply_pfizer.Date,(supply_pfizer.Doses));
grid on;title('Pfizer: supply');
subplot(2,2,3);
plot(supply_jj.Date,cumsum(supply_jj.Doses),'linewidth',1);hold on;
bar(supply_jj.Date,(supply_jj.Doses));
grid on;title('Johnson&Johnson: supply');
subplot(2,2,4);
plot(supply_moderna.Date,cumsum(supply_moderna.Doses),'linewidth',1);hold on;
bar(supply_moderna.Date,(supply_moderna.Doses));
grid on;title('Moderna: supply');

%% vaccination
figure;
plot(vacc_agg.Date,cumsum(vacc_agg.D1),'linewidth',3);hold on;
plot(vacc_agg.Date,cumsum(vacc_agg.D2),'linewidth',3);hold on;
yyaxis right
bar(vacc_agg.Date,(vacc_agg.D1),'facecolor',[0 0 1],'edgecolor',[0 0 1],'facealpha',.25,'edgealpha',.25);hold on;
plot(vacc_agg.Date, movmean(vacc_agg.D1,7),'color',[0 0 .5],'linestyle','-.','linewidth',1);
bar(vacc_agg.Date,(vacc_agg.D2),'facecolor',[1 0 0],'edgecolor',[1 0 1],'facealpha',.25,'edgealpha',.25);hold on;
plot(vacc_agg.Date, movmean(vacc_agg.D2,7),'color',[0.5 0 .5],'linestyle','-.','linewidth',1);
legend({'Dose 1','Dose 2'});
grid on;title('Vaccination progress');

%% prediction
d = 7;
t = 131;
t0 = height(vacc_agg);
waiting_inflow = 21000;
waiting_inflow_decay = [0.01 0.0375 0.065]/d;
loss_d1_d2 = .1;
waiting = 21000.*(1-waiting_inflow_decay.*(0:t)');
waiting_state_d1 = 144835;
waiting_state_d2 = 645569;
offset = ceil(waiting_state_d1/waiting_inflow)-2;
cum_vac_D1 = cumsum(vacc_agg.D1);
cum_vac_D2 = cumsum(vacc_agg.D2);
vaccination_avg_D1 = (cum_vac_D1(end)-cum_vac_D1(end-7))/7;
vaccination_avg_D2 = (cum_vac_D2(end)-cum_vac_D2(end-7))/7;
fcast = table(vacc_agg.Date,vacc_agg.D1,vacc_agg.D1,vacc_agg.D1,vacc_agg.D2,vacc_agg.D2,vacc_agg.D2,vacc_agg.D2*0);
fcast.Properties.VariableNames = {'Date','D1_opt','D1_real','D1_pes','D2_opt','D2_real','D2_pes','supply'};
fcast.Date(end+1:end+t) = fcast.Date(end)+(1:t);
vac_plan_D1 = zeros(t,1);
vac_plan_D1(1:offset) = vaccination_avg_D1;
vac_plan_D1(1+offset:end) = vaccination_avg_D1*(1-waiting_inflow_decay(1).*(0:t-offset-1));
fcast.D1_opt(t0+1:end) = vac_plan_D1;
vac_plan_D1(1+offset:end) = vaccination_avg_D1*(1-waiting_inflow_decay(2).*(0:t-offset-1));
fcast.D1_real(t0+1:end) = vac_plan_D1;
vac_plan_D1(1+offset:end) = max(0,vaccination_avg_D1*(1-waiting_inflow_decay(3).*(0:t-offset-1)));
fcast.D1_pes(t0+1:end) = vac_plan_D1;
vac_plan_D2 = zeros(t,1);
vac_plan_D2(1:offset) = vaccination_avg_D2;
vac_plan_D2(1+offset:end) = vaccination_avg_D2*(1-waiting_inflow_decay(1).*(0:t-offset-1));
fcast.D2_opt(t0+1:end) = vac_plan_D2;
vac_plan_D2(1+offset:end) = vaccination_avg_D2*(1-waiting_inflow_decay(2).*(0:t-offset-1));
fcast.D2_real(t0+1:end) = vac_plan_D2;
vac_plan_D2(1+offset:end) = vaccination_avg_D2*(1-waiting_inflow_decay(3).*(0:t-offset-1));
fcast.D2_pes(t0+1:end) = vac_plan_D2;

opt = 0:0.1:1;opt0=0:.2:1;
pred0 = (opt0.*cumsum(fcast.D1_pes(t0+1:end))+(1-opt0).*cumsum(fcast.D1_real(t0+1:end)));
pred1 = (opt0.*cumsum(fcast.D1_real(t0+1:end))+(1-opt0).*cumsum(fcast.D1_opt(t0+1:end)));
pred=[pred1,pred0(:,2:6)];
idx = find(supply_total.Date(2:end)-supply_total.Date(1:end-1)==0);
idx = idx(end:-1:1);
y = supply_total.Doses;
y(idx)=y(idx)+y(idx+1);y(idx+1) = [];
supply_total(idx+1,:) = [];

supply_total_dates = supply_total.Date;
supply_total_doses = cumsum(supply_total.Doses);
supply_time = days(supply_total.Date-supply_total.Date(1));
tt = t-days(supply_total.Date(end)-today)-1;
tt0 = 20;
supply_total.Doses = cumsum(supply_total.Doses);
supply_val = (supply_total.Doses);
supply_val_0 = interp1(supply_time(end-tt0:end),supply_val((end-tt0:end)),...
    [supply_time(end-tt0:end);(supply_time(end)+1:supply_time(end)+tt)'],'method','extrap');
supply_total.Date(end+1:end+tt) = supply_total.Date(end)+(1:tt);
supply_total.Doses(end-tt+1:end) = supply_val_0(end-tt+1:end);
pred = [zeros(t0,size(pred,2));pred];
cs = cumsum(fcast.D1_opt(1:t0));
pred(1:t0,:) = repmat(cs,1,size(pred,2));
pred(t0+1:end,:) = pred(t0+1:end,:)+cs(end);

figure;
midcolor = .6*[18 181 234]/255+.4*[0 0 1];
T = size(pred,1);
dateFrom = datetime(2020,12,26);
dt = days(fcast.Date(1)-dateFrom);
area(1:t0+dt,supply_total.Doses(end)+zeros(t0+dt,1),...
    'FaceColor',0.75*[1 1 1],'EdgeColor',0.65*[1 1 1],'FaceAlpha',0.5,'EdgeAlpha',0);hold on;
p0=fanChart(1:size(pred,1),pred,pred(:,ceil(size(pred,2)/2)),(opt)','darkcenter',true,'midcolor',midcolor,'alpha',0.35);
p2=plot(1:T,target+0*zeros(T,1),'linestyle','--','linewidth',2,'Color',[.75 .1 0]);hold on;
sup = interp1(days(supply_total.Date-dateFrom+1),supply_total.Doses,1:days(supply_total.Date(end)-dateFrom+1));
plot(sup(dt:T),'k--','linewidth',2);
p1=plot(sup(dt:days(supply_total_dates(end)-dateFrom)),'k','linewidth',2);
p0=plot(cumsum(fcast.D1_opt(1:t0+dt)),'linewidth',2,'Color',midcolor);
ylim([0,supply_total.Doses(end)])
xlim([1,size(pred,1)]);
ylabel('population (mil.)')
legend([p1 p0 p2],{'Vaccines supply','Vaccination: firs dose','First dose vaccination target (approx. 3.3mil)'});
mnths = cumsum([0 31 28 31 30 31 30 31 31 30 31 30 31 31]);
xticks(5+mnths);
xticklabels({'january','february','march','april','may','june','july','august','september','october','november','december','januar'});
title('How many people have had their COVID-19 vaccination?')
grid on;

figure;
plot(fcast.Date,cumsum(fcast.D1_opt),'linewidth',2,'Color',[88 89 91]/255);hold on;
plot(fcast.Date,cumsum(fcast.D1_real),'linewidth',2,'Color',[18 181 234]/255);hold on;
plot(fcast.Date,cumsum(fcast.D1_pes),'linewidth',2,'color',[220 180 123]/255);hold on;

plot(fcast.Date,cumsum(fcast.D2_opt),'linewidth',2,'Color',[88 89 91]/255,'linestyle','-.');hold on;
plot(fcast.Date,cumsum(fcast.D2_real),'linewidth',2,'Color',[18 181 234]/255,'linestyle','-.');hold on;
plot(fcast.Date,cumsum(fcast.D2_pes),'linewidth',2,'color',[220 180 123]/255,'linestyle','-.');hold on;

plot(fcast.Date(1:t0+1),cumsum(fcast.D1_pes(1:t0+1)),'linewidth',2,'color',[.5 .5 .5]);hold on;
plot(supply_total.Date,supply_total.Doses,'k--','linewidth',1);hold on;
plot(supply_total_dates,supply_total_doses,'k','linewidth',1);hold on;
area(datetime(2020,12,26):today,supply_total.Doses(end)+zeros(days(today-datetime(2020,12,26))+1,1),...
    'FaceColor',0.75*[1 1 1],'EdgeColor',0.65*[1 1 1],'FaceAlpha',0.5,'EdgeAlpha',0);hold on;
plot(supply_total.Date,target+0*supply_total.Doses,'linestyle','--','linewidth',2,'Color',[.75 .1 0]);hold on;
xlim([datetime(2020,12,26),fcast.Date(end)])
ylim([0,supply_total.Doses(end)])
ylabel('milions')
grid on;