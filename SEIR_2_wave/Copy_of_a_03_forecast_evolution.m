%% initialization & cleanup
initialize;

%% setup
endFcast = dd(2020,12,31); % end of the year
startWave = dd(2020,09,01); % start of the second wave

%% loading stuff
x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility_new.csv','dateFormat','yyyy-mm-dd','freq','daily');
data = load('inputs.mat','q_mat','Rt','dIt','It','St','s','Rt_last','t0','t1');

%% handle inputs
s = data.s; % s = setparam
startFcast = data.t1+1;
endHist = startFcast-1;
startHist = data.t0;
startEstim = dd(2020,10,10);
startEstimFull = dd(2020,9,15);
delay = 9;
Rt = data.Rt_last;
It = data.It;
St = data.St;

%% handle mobility 
fcastPer = endFcast-startFcast+1;
mobilityFcast = forecast_mobility(mob,fcastPer,startHist,startWave);
Rt_data = tseries(data.t0+1:data.t1,data.q_mat(s.quant_idx_central,:));
mobilityParams = estimate_mobility(mobilityFcast,Rt_data,startEstim,delay,startEstimFull);

%% simulate SEIHR
% simulate augmented SEIR model with clinical section
init.It = It;
init.St = St;
% todo: other init values
% mobility
mobility.values = mobilityFcast.medium(startFcast:end)/100;
mobility.x_grid = mobilityParams.f.x_grid;
mobility.y_grid = mobilityParams.f.y_grid;
mobility.scale = mobilityFcast.medium(startFcast);
% restrictions (NPC as for 15.11.2020 assumed)
restrictions = s.kappa_res_0; % delta, at;
t0 = startFcast;
[res_mean,res_quant] = simulate_SEIHR(fcastPer,Rt,mobility,restrictions,init,t0,s);

%% simulate SIR
alpha_vec = mobilityFcast.medium_norm;
[res_mean,res_quant] = simulate_SIR(fcastPer,Rt,It,St,alpha_vec,endHist,s);

plot_fanchart(res_quant.dIt,s,1,startFcast,endFcast,startFcast);