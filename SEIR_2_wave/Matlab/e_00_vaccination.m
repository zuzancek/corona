%% initialize
initialize;

%% load data
db = dbload('data/vaccination.csv','dateFormat','yyyy-mm-dd','freq','daily');
dbs = dbload('data/supply_plan.csv','dateFormat','yyyy-mm-dd','freq','daily');
today = dd(2021,3,19-1);
dateFrom = startdate(db.Cum);
dateTo_short = dd(2021,4,1);
dateTo = dd(2021,12,31);
target = 3300000;

%% process inputs
supply_short = db.Cum;
supply_realised = resize(db.Cum,dateFrom:today);
supply_plan_short = resize(db.Cum,today:dateTo_short);
supply_short_tnd = hpf(supply_short,'lambda',14400);
supply_short_tnd(find(supply_short_tnd<0)) = 0;
d1_cum = cumsum(db.D1);
d2_cum = cumsum(db.D2);
covered = d2_cum;
covered_realised = resize(covered,dateFrom:today);
d_cum = d1_cum+d2_cum;
vaccination_realised = resize(d_cum,dateFrom:today);

%% short-term analysis
% ******* Vaccination in Q1
d1_cum_tnd = hpf(resize(d1_cum,dateFrom:today),'lambda',14400);
d1_cum_tnd(find(d1_cum_tnd<0)) = 0; %#ok<*FNDSB>
d2_cum_tnd = hpf(resize(d2_cum,dateFrom:today),'lambda',14400);
d2_cum_tnd(find(d2_cum_tnd<0)) = 0;
d_cum_tnd = hpf(resize(d_cum,dateFrom:today),'lambda',14400);
d_cum_tnd(find(d_cum_tnd<0)) = 0;

% plot results
figure;
subplot(2,1,1)
bar(supply_short,'FaceAlpha',0.8); hold on;
bar(resize(d_cum,dateFrom:today),'FaceAlpha',0.8)
plot(supply_short_tnd,'b','linewidth',2,'linestyle','--');
plot(d_cum_tnd,'r','linewidth',2,'linestyle','--');
grid on;
xlim([dateFrom dateTo_short]);
ylim([0 round(supply_short(end),2,'significant')])
legend({'dodane', 'pouzite'})
title('Dodavky vakcin v 2021Q1');

subplot(2,1,2)
bar(d_cum,'FaceColor',0.5*[1 1 1],'EdgeColor',0.75*[1 1 1],'FaceAlpha',0.8); hold on;
plot(d1_cum,'b','linewidth',2); hold on;
plot(d2_cum,'r','linewidth',2); hold on;
plot(d_cum_tnd,'k','linewidth',2,'linestyle','-');
grid on;
xlim([dateFrom today]);
legend({'celkovo podanych davok','prva davka','druha davka (=plna ochrana)'})
title('Priebeh vakcinacie v 2021Q1');

%% long-term analysis
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

% ********* supply plan 
dbs.astra(today) = sum(double(db.astra),'omitnan'); dbs.astra = interp(dbs.astra);
dbs.pfizer(today) = sum(double(db.pfizer),'omitnan');   dbs.pfizer = interp(dbs.pfizer); 
dbs.moderna(today) = sum(double(db.moderna),'omitnan'); dbs.moderna = interp(dbs.moderna);
dbs.johnson(today) = 0;      dbs.johnson = interp(dbs.johnson);
dbs.curevac(today) = 0;      dbs.curevac = interp(dbs.curevac);
supply_plan = dbs.astra+dbs.moderna+dbs.pfizer+dbs.curevac+dbs.johnson;
supply_plan_cover = 0.5*(dbs.astra+dbs.moderna+dbs.pfizer+dbs.curevac)+dbs.johnson;
supply_plan(dateFrom+1:dateTo_short-5) = NaN; supply_plan = interp(supply_plan);
supply_plan_cover(dateFrom+1:dateTo_short-5) = NaN; supply_plan_cover = interp(supply_plan_cover);
supply_plan_tnd = hpf(supply_plan,'lambda',14400);

% ********** vaccination plan (multiple scenaria)
d_vac = vaccination_realised(today)/supply_realised(today);
vaccination_plan_0 = d_vac*supply_plan; vaccination_plan_0(dateFrom:today) = d_cum_tnd;
vaccination_plan_1 = supply_plan;       vaccination_plan_1(dateFrom:today) = d_cum_tnd;
vaccination_plan = 0.5*(vaccination_plan_0+vaccination_plan_1);
opt = 0:0.1:1;
vaccination_plan_opt = (opt.*vaccination_plan_0+(1-opt).*vaccination_plan_1)';
vaccination_plan_opt(today+1:today+21,:) = NaN;
vaccination_plan_opt = interp(vaccination_plan_opt);
vp = double(vaccination_plan_opt);
vaccination_plan(today+1:today+21) = NaN;
vaccination_plan = interp(vaccination_plan);
d_cov = covered_realised(today)*2.75/supply_realised(today);
cover_plan_0 = d_cov*supply_plan_cover; cover_plan_0(dateFrom:today) = d2_cum_tnd;
cover_plan_1 = supply_plan_cover;       cover_plan_1(dateFrom:today) = d2_cum_tnd;
cover_plan = 0.5*(cover_plan_0+cover_plan_1);
cover_plan_opt = (opt.*cover_plan_0+(1-opt).*cover_plan_1)';
cover_plan_opt(today+1:today+21,:) = NaN;
cover_plan_opt = interp(cover_plan_opt);
cp = double(cover_plan_opt);
cover_plan(today+1:today+21) = NaN;
cover_plan = interp(cover_plan);

% plot results
figure;
hist = total_num+zeros(today-dateFrom+1,1);
pp1=plot(double(supply_plan),'linewidth',2); hold on;
pp2=plot(double(vaccination_plan),'linewidth',2); hold on;
fanChart(1:size(vp,1),vp,vp(:,ceil(size(vp,2)/2)),(opt)','darkcenter',true,'midcolor',pp2.Color,'alpha',0.75);
area(hist,'FaceColor',0.75*[1 1 1],'EdgeColor',0.65*[1 1 1],'FaceAlpha',0.5,'EdgeAlpha',0);hold on;
plot(double(supply_plan),'linewidth',2,'color',pp1.Color); hold on;
pp3=plot(double(supply_realised),'linewidth',2,'Color',0.5*[0 0.5 1]);
pp4=plot(double(vaccination_realised),'linewidth',2, 'Color',0.5*[1 0.5 0]);
grid on;
xlim([0 length(supply_plan)-2])
mnths = cumsum([0 31 28 31 30 31 30 31 31 30 31 30 31 31]);
xticks(5+mnths);
xticklabels({'januar','februar','marec','april','maj','jun','jul','august','september','oktober','november','december','januar'});
ylim([0 total_num]);
ylabel('kusy vakcin (mil.)');
legend([pp1 pp2 pp3 pp4],{'plan dodavky vakcin','planovana realizacia vakcinacie','realizacia dodavok','realizacia vakcinacie'});
title('Priebeh vakcinacie v roku 2021 : Pocty podanych vakcin');

figure;
hist = total_num+zeros(today-dateFrom+1,1);
pp1=plot(double(supply_plan_cover),'linewidth',2); hold on;
pp2=plot(double(cover_plan),'linewidth',2); hold on;
fanChart(1:size(cp,1),cp,cp(:,ceil(size(cp,2)/2)),(opt)','darkcenter',true,'midcolor',pp2.Color,'alpha',0.75);
area(hist,'FaceColor',0.75*[1 1 1],'EdgeColor',0.65*[1 1 1],'FaceAlpha',0.5,'EdgeAlpha',0);hold on;
plot(double(supply_plan_cover),'linewidth',2,'color',pp1.Color); hold on;
plot(double(resize(supply_realised,dateFrom:today)/2),'linewidth',2,'Color',0.5*[0 0.5 1]);
plot(double(resize(d2_cum,dateFrom:today)),'linewidth',2, 'Color',0.5*[1 0.5 0]);
pp3=plot(target+zeros(dateTo-dateFrom+1,1),'k--','linewidth',2);
grid on;
ylabel('populacia (mil.)')
xlim([0 length(supply_plan)-2])
mnths = cumsum([0 31 28 31 30 31 30 31 31 30 31 30 31 31]);
xticks(5+mnths);
xticklabels({'januar','februar','marec','april','maj','jun','jul','august','september','oktober','november','december','januar'});
ylim([0 7*10^6]);
legend([pp1 pp2 pp3],{'plan dodavok (1+2 doza)','planovana zaockovanost populacie','ciel zaockovanosti (60% populacie)'});
title('Dosiahnutie ciela v zaockovanosti (60% populacie, cca 3.3mil), rok 2021');

figure;
cum_supply = zeros(5,length(dbs.astra));
cum_supply(1,:) = dbs.pfizer;
cum_supply(2,:) = dbs.moderna;
cum_supply(3,:) = dbs.astra;
cum_supply(4,:) = dbs.curevac;
cum_supply(5,:) = dbs.johnson;
bar(cum_supply','stacked','FaceAlpha',1);
grid on;
xlim([0 length(supply_plan)-2]);
ylabel('kusy vakcin (mil.)');
mnths = cumsum([0 31 28 31 30 31 30 31 31 30 31 30 31]);
xticks(5+mnths);
xticklabels({'januar','februar','marec','april','maj','jun','jul','august','september','oktober','november','december'});
legend({'BioNTech/Pfizer','Moderna','Oxford/AstraZeneca','Curevac','Johnson&Johnson'});
title('Plan dodavok vakcin v roku 2021');
