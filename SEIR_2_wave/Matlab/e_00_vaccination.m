initialize;
today = datetime(2021,05,24);

%% load
supply_astra=readtable('data/supply_plan_astrazeneca.csv');
supply_pfizer=readtable('data/supply_plan_pfizer.csv');
supply_jj=readtable('data/supply_plan_jj.csv');
supply_moderna=readtable('data/supply_plan_moderna.csv');
supply_total = sortrows(union(union(union(supply_astra,supply_pfizer,'stable'),supply_jj,'stable'),supply_moderna,'stable'),1);
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
