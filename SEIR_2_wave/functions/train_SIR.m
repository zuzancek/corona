function [res_mean,res_quant] = train_SIR(time_interval,assumptions,inputs,s)

%% handle inputs
dateFrom = time_interval.dateFrom;
dateTo = time_interval.dateTo;
T = dateTo-dateFrom+1;
% assumptions - paths
mobility = assumptions.mobility;
restrictions = assumptions.restrictions;
alpha = assumptions.asymp_ratio;
tau = assumptions.obs_ratio;
obs_ratio_effect = assumptions.obs_ratio_effect;
% inputs: Rt, init values
init = inputs.init;
Rt = inputs.Rt;
% other parameters (stored in s)
q_vec = s.quant;
pop_size = s.pop_size;
N = s.sim_num;
% alpha_hr = s.alpha_hr;
lambda = s.lambda;
alpha_hd = 9.71/100; %s.alpha_hd;
self_isolation_effect = s.self_isolation_effect;
case_isolation_effect = s.case_isolation_effect;

% init values (observed)
St = init.St(1);
Iot = init.Iot(1);
if isfield(init,'Iat'); Iat = init.Iat(1); else; Iat = alpha(1)*Iot; end
if isfield(init,'Ist'); Ist = init.Ist(1); else; Ist = (1-alpha(1))*Iot; end
Ht = init.Ht(1);
Dt = init.Dt(1);

% setup
S_vec = zeros(N,T+1);       S_vec(:,1) = St;
Iu_vec = zeros(N,T+1);      Iu_vec(:,1) = (1/tau(1)-1)*Iot;
Io_vec = zeros(N,T+1);      Io_vec(:,1) = Iot;
Ia_vec = zeros(N,T+1);      Ia_vec(:,1) = Iat; 
Is_vec = zeros(N,T+1);      Is_vec(:,1) = Ist; 
H_vec = zeros(N,T+1);       H_vec(:,1) = Ht;
D_vec = zeros(N,T+1);       D_vec(:,1) = Dt;
R_vec = zeros(N,T+1);
dI_in_vec = zeros(N,T+1);     
R_eff = zeros(N,T+1);       R_eff(:,1) = Rt(:,1);
st = zeros(T+1,1);
it = st; iot = st; iut = st; ist = st; iat = st; dt = st; rt = st; ht = st;

%% initialize
T_si_vec = get_rv(s.SI);
T_hosp_vec = get_rv(s.T_hosp);
T_death.mean = s.T_death.mean; T_death.std = 1;
T_death = get_rv(T_death);
T_rec.mean = s.T_rec; T_rec.std = 1;
T_rec = get_rv(T_rec);
shift = min(floor(get_rv(s.SI)),ceil(2.5*s.SI.mean));
shift_vec = [1:N]'-shift.*N; %#ok<NBRAK>
gamma = 1./T_si_vec;
gamma_hosp_vec = 1./T_hosp_vec;

idx = ones(N,1);
kappa_res = get_kappa_res();
kappa_mob = get_kappa_mob();

varsigma_unobs = 1/self_isolation_effect*ones(T,1);
varsigma = ((s.T_inf_obs.mean-s.T_inf_obs0.mean)+s.T_inf_obs0.mean/case_isolation_effect)/s.T_inf_unobs.mean;

for t=1:T
    r_idx = (t+dateFrom-1).*N+shift_vec;
    R_eff(:,t) = Rt(r_idx).*kappa_mob(t).*kappa_res(t).*obs_ratio_effect(t);
    dI_in_vec(:,t) = R_eff(:,t).*S_vec(:,t)/pop_size.*...
        (varsigma.*Io_vec(:,t)+varsigma_unobs(t).*Iu_vec(:,t)).*gamma;
    S_vec(:,t+1) = S_vec(:,t)-dI_in_vec(:,t);   
    Io_vec(:,t) = Ia_vec(:,t)+Is_vec(:,t);
    Iu_vec(:,t+1) = Iu_vec(:,t)+(1-s.obs_ratio).*dI_in_vec(:,t)-gamma.*Iu_vec(:,t);
    Ia_vec(:,t+1) = Ia_vec(:,t)+alpha(t).*s.obs_ratio.*dI_in_vec(:,t)-gamma.*Ia_vec(:,t); 
    Is_vec(:,t) = Is_vec(:,t)+(1-alpha(t)).*s.obs_ratio.*dI_in_vec(:,t)...
        -(1-lambda).*gamma.*Is_vec(:,t)-lambda.*gamma_hosp_vec.*Is_vec(:,t); 
    H_vec(:,t+1) = H_vec(:,t).*(1-(1-alpha_hd)./T_rec-alpha_hd./T_death)+lambda.*gamma_hosp_vec.*Is_vec(:,t);
    D_vec(:,t+1) = D_vec(:,t)+alpha_hd./T_death.*H_vec(:,t);
    R_vec(:,t+1) = R_vec(:,t)+(1-alpha_hd)./T_rec.*H_vec(:,t)+gamma.*Ia_vec(:,t)...
        +(gamma-lambda.*gamma_hosp_vec./(1-alpha(t))).*Is_vec(:,t);
    idx = idx & Io_vec(:,t+1)>=0 & Iu_vec(:,t+1)>=0 & H_vec(:,t+1)>=0;
end

idx = find(idx>0);
S_vec = S_vec(idx,:);
Iu_vec = Iu_vec(idx,:);
Io_vec = Io_vec(idx,:);
Ia_vec = Ia_vec(idx,:);
Is_vec = Is_vec(idx,:);
H_vec = H_vec(idx,:);
D_vec = D_vec(idx,:);
R_vec = R_vec(idx,:);
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
    ht(t) = mean(H_vec(:,t));
    dt(t) = mean(D_vec(:,t));
    rt(t) = mean(R_vec(:,t));
end
res_mean.It = it;
res_mean.Iat = iat;
res_mean.Ist = ist;
res_mean.Iot = iot;
res_mean.Iut = iut;
res_mean.St = st;
res_mean.Rt = rt;
res_mean.Dt = dt;
res_mean.Ht = ht;

M = length(q_vec);
res_quant.Iat = get_quant(Ia_vec(:,1:T));
res_quant.Ist = get_quant(Is_vec(:,1:T));
res_quant.Iut = get_quant(Iu_vec(:,1:T));
res_quant.Iot = get_quant(Io_vec(:,1:T));
res_quant.It = get_quant(I_vec(:,1:T));
res_quant.St = get_quant(S_vec(:,1:T));
res_quant.Dt = get_quant(D_vec(:,1:T));
res_quant.Ht = get_quant(H_vec(:,1:T));
res_quant.Rt = get_quant(R_vec(:,1:T));

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