%% initialization & cleanup
initialize;

%% setup
endFcast = dd(2020,12,31); % end of the year
startWave = dd(2020,09,01); % start of the second wave

%% loading stuff
x = dbload('data/korona_data.csv','dateFormat','yyyy-mm-dd','freq','daily');
mob = dbload('data/mobility.csv','dateFormat','yyyy-mm-dd','freq','daily');
data = load('results_Rt_SIR_pcr.mat','s','t0','t1','q_mat','Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
data_real = load('results_Rt_SIR_real_pcr.mat','s','t0','t1','q_mat','Rt','Yt','Rt_last','Rt_dist','Rt_rnd');
mf = load('mobility_forecast.mat','mobilityFcast','startHist','fcastPer','startWave');
me = load('mobility_estimation.mat','mobilityParams','mobilityParams_real','startEstim','delay','startEstimFull');

%% handle inputs
% ****************** Dates
% ****************** Mobility
% medium scenario assumed
mobility.forecast = true;
mobility.values = mf.mobilityFcast.medium(startFcast:end)/100;
mobility.scale = mf.mobilityFcast.medium(startFcast);
mobility.x_grid = me.mobilityParams.x;
mobility.y_grid = me.mobilityParams.y_pos;


s = data.s; % s = setparam
startFcast = data.t1+1;
endHist = startFcast-1;
startHist = data.t0;
startEstim = dd(2020,10,10);
startEstimFull = dd(2020,9,15);
endFcast = dd(2020,12,31);

%% simulate SEIHR
% simulate augmented SEIR model with clinical section
init = data.Yt;
% mobility
mobility.values = mf.mobilityFcast.medium(startFcast:end)/100;
mobility.x_grid = me.mobilityParams.f.x_grid;
mobility.y_grid = me.mobilityParams.f.y_grid;
mobility.scale = mf.mobilityFcast.medium(startFcast);
% restrictions (NPC as for 15.11.2020 assumed)
restrictions = s.kappa_res_0; % delta, at;
t0 = startFcast;
Rt = data.Rt_last;
fcastPer = endFcast-startFcast+1;
[res_mean,res_quant] = simulate_SEIHR(fcastPer,Rt,mobility,restrictions,init,t0,s);

%% plot results
plot_fanchart(res_quant.dIt,s,1,startFcast,endFcast,startFcast);