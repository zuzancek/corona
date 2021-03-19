initialize;

db = dbload('data/vaccination.csv','dateFormat','yyyy-mm-dd','freq','daily');
dbs = dbload('data/supply_plan.csv','dateFormat','yyyy-mm-dd','freq','daily');
today = dd(2021,3,19);
dateFrom = startdate(db.Cum);
dateTo_short = dd(2021,4,1);
dateTo = dd(2021,12,31);

supply_short = db.Cum;
supply_realised = resize(db.Cum,dateFrom:today);
supply_plan_short = resize(db.Cum,today:dateTo_short);
supply_short_tnd = hpf(supply_short,'lambda',14400);
supply_short_tnd(find(supply_short_tnd<0)) = 0;
d1_cum = cumsum(db.D1);
d2_cum = cumsum(db.D2);
d_cum = d1_cum+d2_cum;
vaccination_realised = resize(d_cum,dateFrom:today);

% ******* Vaccination in Q1
d1_cum_tnd = hpf(resize(d1_cum,dateFrom:today),'lambda',14400);
d1_cum_tnd(find(d1_cum_tnd<0)) = 0; %#ok<*FNDSB>
d2_cum_tnd = hpf(resize(d2_cum,dateFrom:today),'lambda',14400);
d2_cum_tnd(find(d2_cum_tnd<0)) = 0;
d_cum_tnd = hpf(resize(d_cum,dateFrom:today),'lambda',14400);
d_cum_tnd(find(d_cum_tnd<0)) = 0;

figure;
subplot(2,1,1)
bar(supply_short); hold on;
bar(resize(d_cum,dateFrom:today))
plot(supply_short_tnd,'b','linewidth',1,'linestyle','--');
plot(d_cum_tnd,'r','linewidth',1,'linestyle','--');
grid on;
xlim([dateFrom dateTo_short]);
ylim([0 round(supply_short(end),2,'significant')])
legend({'dodane', 'pouzite'})
title('Dodavky vakcin v 2021Q1');

subplot(2,1,2)
bar(d_cum,'FaceColor',0.5*[1 1 1],'EdgeColor',0.75*[1 1 1]); hold on;
plot(d1_cum,'b','linewidth',1); hold on;
plot(d2_cum,'r','linewidth',1); hold on;
plot(d_cum_tnd,'k','linewidth',2,'linestyle','-');
grid on;
xlim([dateFrom today]);
legend({'celkovo podanych davok','prva davka','druha davka'})
title('Priebeh vakcinacie v 2021Q1');

% **** contract plan
plan.pfizer.contract = 2407086;
plan.pfizer.mult = 1;
plan.astrazeneca.contract = 3638830;
plan.astrazeneca.mult = 1;
plan.moderna.contract = 962917;
plan.moderna.mult = 1;
plan.johnson.contract = 2407086;
plan.johnson.mult = 2;
total_num = plan.pfizer.contract.*plan.pfizer.mult+plan.astrazeneca.contract*plan.astrazeneca.mult+...
    plan.moderna.contract*plan.moderna.mult+plan.johnson.contract*plan.johnson.mult;

% % supply plan - worst case scen.
% supply_plan_wcs = tseries(dateFrom:dateTo,NaN);
% supply_plan_wcs(dateFrom:dateTo_short) = supply_short_tnd;
% supply_plan_wcs(dateTo) = total_num;
% supply_plan_wcs = interp(supply_plan_wcs);
% x0=today-dateFrom; y0 = supply_plan_wcs(today);
% x1=dateTo-dateFrom; y1 = supply_plan_wcs(dateTo);
% q = log(y1/y0)/log(x1/x0); p = y0/x0^q;
% supply_plan_an = p*(0:x1)'.^q+0*supply_plan_wcs;

% supply plan 
dbs.astra(today) = sum(double(db.astra),'omitnan');     
plan.astrazeneca.data = interp(dbs.astra);
dbs.pfizer(today) = sum(double(db.pfizer),'omitnan');   plan.pfizer.data = interp(dbs.pfizer);
dbs.moderna(today) = sum(double(db.moderna),'omitnan'); plan.moderna.data = interp(dbs.moderna);
dbs.johnson(today) = 0;      plan.johnson.data = interp(dbs.johnson);
dbs.curevac(today) = 0; plan.curevac.data = interp(dbs.curevac);
supply_plan = plan.astrazeneca.data+plan.moderna.data+plan.pfizer.data+plan.curevac.data+plan.johnson.data;
supply_plan(dateFrom:today) = resize(supply_short,dateFrom:today);
supply_plan_tnd = hpf(supply_plan,'lambda',14400);

vaccination_plan = tseries(dateFrom:dateTo,NaN);
vaccination_plan(dateFrom:today) = d_cum_tnd;
vaccination_plan(dateTo) = 0;
x0=today-dateFrom-14; y0 = vaccination_plan(today-14);
x1=today-dateFrom; y1 = vaccination_plan(today);
q = log(y1/y0)/log(x1/x0); p = y0/x0^q;
vaccination_plan_an = p*(0:dateTo-dateFrom)'.^q+0*vaccination_plan;
vaccination_plan(dateTo) = vaccination_plan_an;
vaccination_plan = interp(vaccination_plan);

figure;
hist = total_num+tseries(dateFrom:today,0);
area(hist,'FaceColor',0.75*[1 1 1],'EdgeColor',0.75*[1 1 1],'FaceAlpha',0.25,'EdgeAlpha',0);hold on;
pp1=plot(supply_plan,'linewidth',2,'Color',[0 0.5 1]); hold on;
pp2=plot(vaccination_plan,'linewidth',2,'Color',[1 0.5 0]); hold on;
pp3=plot(supply_realised,'linewidth',2,'Color',0.5*[0 0.5 1]);
pp4=plot(vaccination_realised,'linewidth',2, 'Color',0.5*[1 0.5 0]);
% plot(supply_plan_an,'linewidth',1,'Color',0.75*[1 1 1]); hold on;
% plot(vaccination_plan_an,'linewidth',1,'Color',0.5*[1 1 1]); hold on;
grid on;
xlim([dateFrom dateTo]);
ylim([0 total_num]);
legend([pp1 pp2 pp3 pp4],{'plan dodavky vakcin','planovana realizacia vakcinacie','realizacia dodavok','realizacia vakcinacie'});
title('Priebeh vakcinacie v roku 2021');

