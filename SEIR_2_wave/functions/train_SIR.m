function [res_mean,res_quant] = train_SIR(time_interval,inputs,s)

%% handle inputs
dateFrom = time_interval.dateFrom;
dateTo = time_interval.dateTo;
T = dateTo-dateFrom+1;
mobility = inputs.mobility;
restrictions = inputs.restrictions;
q_vec = s.quant;
obs_ratio_effect = inputs.obs_ratio_effect;
pop_size = s.pop_size;
init = inputs.init;
Rt = inputs.Rt;
N = s.sim_num;
% T_test0 = inputs.T_test;
alpha = inputs.asymp_ratio;
tau = inputs.obs_ratio;
% lambda = s.lambda;
s.share_reas =1;

% init values (observed)
St = init.S(1);
Iot = init.Io(1);

% setup
S_vec = zeros(N,T+1);       S_vec(:,1) = St;
Iu_vec = zeros(N,T+1);      Iu_vec(:,1) = (1/tau(1)-1)*Iot;
Io_vec = zeros(N,T+1);      Io_vec(:,1) = Iot;
Ia_vec = zeros(N,T+1);      Ia_vec(:,1) = alpha(1)*Iot; 
Is_vec = zeros(N,T+1);      Is_vec(:,1) = (1-alpha(1))*Iot; 
dI_in_vec = zeros(N,T+1);     
R_eff = zeros(N,T+1);       R_eff(:,1) = Rt(:,1);
st = zeros(T+1,1);
it = st; iot = st; iut = st; ist = st; iat = st;

%% initialize
T_si_vec = get_rv(s.SI);
shift = min(floor(get_rv(s.SI)),ceil(2.5*s.SI.mean));
shift_vec = [1:N]'-shift.*N;
gamma = 1./T_si_vec;

% T_pre = get_rv(s.T_pre);
% T_test = repmat(T_test0',N,1)+repmat(T_pre,1,T);
% T_inf_obs0 = get_rv(s.T_inf_obs0); 
% T_inf_obs = T_inf_obs0+T_test;
% gamma_lat = 1./T_lat;
% gamma_unobs = 1./T_inf_unobs;
% gamma_obs = 1./T_inf_obs;
% gamma_obs0 = 1./T_inf_obs0;
% gamma_test = 1./T_test;
% gamma_hosp = 1./T_hosp;

idx = ones(N,1);
kappa_res = get_kappa_res();
kappa_mob = get_kappa_mob();

s.self_isolation_effect = 1-0.05;
varsigma_unobs = 1/s.self_isolation_effect*ones(T,1);

varsigma = ((s.T_inf_obs.mean-s.T_inf_obs0.mean)+s.T_inf_obs0.mean/s.case_isolation_effect)/s.T_inf_unobs.mean;

for t=1:T
    r_idx = (t+dateFrom-1).*N+shift_vec;
    R_eff(:,t) = Rt(r_idx).*kappa_mob(t).*kappa_res(t).*obs_ratio_effect(t);
    dI_in_vec(:,t) = R_eff(:,t).*S_vec(:,t)/pop_size.*...
        (varsigma.*Io_vec(:,t)+varsigma_unobs(t).*Iu_vec(:,t)).*gamma;
    S_vec(:,t+1) = S_vec(:,t)-dI_in_vec(:,t);   
    Io_vec(:,t+1) = Io_vec(:,t)+s.obs_ratio.*dI_in_vec(:,t)-gamma.*Io_vec(:,t);
    Iu_vec(:,t+1) = Iu_vec(:,t)+(1-s.obs_ratio).*dI_in_vec(:,t)-gamma.*Iu_vec(:,t);
    Ia_vec(:,t) = alpha(t).*Io_vec(:,t); Is_vec(:,t) = Io_vec(:,t)-Ia_vec(:,t);
    idx = idx & Io_vec(:,t+1)>=0 & Iu_vec(:,t+1)>=0;
end

idx = find(idx>0);
S_vec = S_vec(idx,:);
Iu_vec = Iu_vec(idx,:);
Io_vec = Io_vec(idx,:);
Ia_vec = Ia_vec(idx,:);
Is_vec = Is_vec(idx,:);
I_vec = Io_vec+Iu_vec;

% mean values
res_mean = struct;
for t = 1:T+1
    it(t) = mean(I_vec(:,t));
    st(t) = mean(S_vec(:,t));
    iat(t) = mean(Ia_vec(:,t));
    ist(t) = mean(Is_vec(:,t));
    iot(t) = mean(Io_vec(:,t));
    iut(t) = mean(Iu_vec(:,t));
end
res_mean.It = it;
res_mean.Iat = iat;
res_mean.Ist = ist;
res_mean.Iot = iot;
res_mean.Iut = iut;
res_mean.St = st;

M = length(q_vec);
res_quant.Iat = get_quant(Ia_vec(:,1:T));
res_quant.Ist = get_quant(Is_vec(:,1:T));
res_quant.Iut = get_quant(Iu_vec(:,1:T));
res_quant.Iot = get_quant(Io_vec(:,1:T));
res_quant.It = get_quant(I_vec(:,1:T));
res_quant.St = get_quant(S_vec(:,1:T));

figure('Name','Observed');
plot(init.Io);hold on;
plot(res_mean.Iot);hold on;
plot(res_quant.Iot(10,:));hold on;
grid on;
legend({'data implied','sim mean','sim median'});

% figure;
% plot(Rt);hold on;
% plot(R_eff);
% grid on;

%% helpers

    function [x] = get_rv(y)
        shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
        L = length(shape0);
        shape0_vec = repmat(shape0,N,1);
        scale0_vec = scale0*ones(N,L);
        x = gamrnd(shape0_vec,scale0_vec);
    end

    function [x] = get_quant(vec)
        x = zeros(M,T);
        for j = 1:M
            x(j,:) = quantile(vec,q_vec(j),1);
        end
    end

    function [r] = get_kappa_res()
        if restrictions.forecast
            dr = s.kappa_res_delta_0+restrictions.delta*min([0:T],restrictions.at)/restrictions.at; %#ok<NBRAK>
            r = (1+dr.^s.kappa_res_alpha).*s.beta_res;
        else
            r = ones(T,1);
        end
    end

    function [m] = get_kappa_mob()
        if mobility.forecast
            m0 = interp1(mobility.x_grid,mobility.y_grid,mobility.values);
            ms = interp1(mobility.x_grid,mobility.y_grid,mobility.scale);
            m = m0/ms;
        else
            m = ones(T,1);
        end
    end

    

end