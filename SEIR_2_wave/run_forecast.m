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
Rt = data.Rt_last;
It = data.It;
St = data.St;

%% forecast mobility evolution
fcastPer = endFcast-startFcast+1;
mobilityFcast = estimate_mobility(mob,fcastPer,startHist,startFcast-1,startWave);

%% simulate SIR
alpha_vec = mobilityFcast.low_norm;
[res_mean,res_quant] = simulate_SIR(fcastPer,Rt,It,St,alpha_vec,endHist,s);

plot_fanchart(res_quant.dIt,s,1,startFcast,endFcast,startFcast);