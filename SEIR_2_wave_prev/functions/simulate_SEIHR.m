function [res_mean,res_q] = simulate_SEIHR(T,Rt,mobility,restrictions,init,t0,s)

% init values (observed)
Iat = init.Ia;
Ist = init.Is;
Iut = init.Iu;
Et = init.E;
St = init.S;
Ht = init.H;
Dt = init.D;
Vt = init.V;

% initialize
T_inf_asymp = s.T_inf_asymp;
T_inf_symp = s.T_inf_symp;
T_inf_unobs = s.T_inf_unobs;
T_inf_hosp = s.T_inf_hosp_vec;
T_lat = s.T_lat;
T_pre = s.T_pre;
rho = s.obs_ratio;
% sigma = s.symp_ratio_obs;
theta = s.p_a_s;
lambda = s.lambda;
omega_novent = s.omega_novent; 
omega_vent = s.omega_vent;
psi_novent = s.psi_novent; 
psi_vent = s.psi_vent;
xi = s.xi;
iota = s.iota;
pop_size = s.pop_size;
gamma_novent = s.gamma_novent;
varsigma = s.self_isolation_effect;
q_vec = s.q_vec;

% setup
N = length(Rt);
T_inf_asymp_vec = get_rv(T_inf_asymp);
T_inf_symp_vec = get_rv(T_inf_symp);
T_inf_unobs_vec = get_rv(T_inf_unobs);
T_inf_hosp_vec = get_rv(T_inf_hosp);
T_lat_vec = get_rv(T_lat);
T_pre_vec = get_rv(T_pre);

gamma_lat = 1./T_lat_vec;
gamma_pre = 1./T_pre_vec;
gamma_asymp = 1./T_inf_asymp_vec;
gamma_symp = 1./T_inf_symp_vec;
gamma_unobs = 1./T_inf_unobs_vec;
gamma_hosp = 1./T_inf_hosp_vec;

S_vec = zeros(N,T+1);       S_vec(:,1) = St(t0);
E_vec = zeros(N,T+1);       E_vec(:,1) = Et(t0);
Ia_vec = zeros(N,T+1);      Ia_vec(:,1) = Iat(t0); % rho*(1-sigma)*It(t0);
Is_vec = zeros(N,T+1);      Is_vec(:,1) = Ist(t0); % rho*sigma*It(t0);
Iu_vec = zeros(N,T+1);      Iu_vec(:,1) = Iut(t0); % (1-rho)*It(t0);
R_vec = zeros(N,T+1);
H_vec = zeros(N,T+1);       H_vec(:,1) = Ht(t0);
V_vec = zeros(N,T+1);       V_vec(:,1) = Vt(t0);
N_vec = zeros(N,T+1);       N_vec(:,1) = Ht(t0)-Vt(t0);
D_vec = zeros(N,T+1);       D_vec(:,1) = Dt(t0);
dI_in_vec = zeros(N,T+1);   dI_in_vec(:,1) = Rt;
dE_in_vec = zeros(N,T+1);   dE_in_vec(:,1) = Rt;
R_eff = zeros(N,T+1);       R_eff(:,1) = Rt;

idx = ones(N,1);
kappa_res = get_kappa_res();
kappa_mob = get_kappa_mob();
% calculation
for t=1:T
    % epidemiological part
    R_eff(:,t+1) = Rt.*kappa_mob(t).*kappa_res(t);
    dE_in_vec(:,t) = R_eff(:,t).*S_vec(:,t)/pop_size.*...
        (Is_vec(:,t).*((1-lambda)*gamma_symp+lambda*gamma_hosp)+...
        Ia_vec(:,t).*gamma_asymp+...
        Iu_vec(:,t).*gamma_unobs/varsigma);
    S_vec(:,t+1) = S_vec(:,t)-dE_in_vec(:,t);
    E_vec(:,t+1) = E_vec(:,t).*(1-gamma_lat)+dE_in_vec(:,t);
    Iu_vec(:,t+1) = Iu_vec(:,t)+(1-rho).*gamma_lat.*E_vec(:,t)-gamma_unobs.*Iu_vec(:,t);
    Ia_vec(:,t+1) = Ia_vec(:,t)+rho.*gamma_lat.*E_vec(:,t)-theta.*gamma_pre.*Ia_vec(:,t)...
        -(1-theta).*gamma_asymp.*Ia_vec(:,t);    
    Is_vec(:,t+1) = Is_vec(:,t)+theta.*gamma_pre.*Ia_vec(:,t)-lambda*gamma_hosp.*Is_vec(:,t)...
        -(1-lambda).*gamma_symp.*Is_vec(:,t);
    % clinical part
    N_vec(:,t+1) = lambda*gamma_hosp.*Is_vec(:,t)+N_vec(:,t).*(omega_novent.*(1-psi_novent)...
        +(1-omega_novent).*(1-gamma_novent-xi));
    V_vec(:,t+1) = (1-omega_novent).*xi.*N_vec(:,t)+V_vec(:,t).*(omega_vent.*(1-psi_vent)+(1-omega_vent).*(1-gamma_vent));
    H_vec(:,t+1) = N_vec(:,t+1)+V_vec(:,t+1);
    % final stages
    R_vec(:,t+1) = R_vec(:,t)+gamma_unobs.*Iu_vec(:,t)+(1-theta).*gamma_asymp.*Ia_vec(:,t)...
        +(1-lambda).*gamma_symp.*Is_vec(:,t)+(1-omega_novent).*gamma_novent.*N_vec(:,t)...
        +(1-omega_vent).*gamma_vent.*V_vec(:,t);
    D_vec(:,t+1) = D_vec(:,t)+omega_novent.*psi_novent.*N_vec(:,t)+omega_vent.*psi_vent.*V_vec(:,t);
    % correct indices only
    idx = idx & Ia_vec(:,t+1)>0 & Is_vec(:,t+1)>0 & Iu_vec(:,t+1)>0;
end

% adjust for valid indices
idx = find(idx>0);
S_vec = S_vec(idx,:);
Ia_vec = Ia_vec(idx,:);
Is_vec = Is_vec(idx,:);
Iu_vec = Iu_vec(idx,:);
Io_vec = Ia_vec+Is_vec;
dI_in_vec = dI_in_vec(idx,:);
N_vec = N_vec(idx,:);
C_vec = iota*N_vec;
N_vec = N_vec-C_vec;
C_vec = C_vec(idx,:);
V_vec = V_vec(idx,:);
H_vec = H_vec(idx,:);
R_vec = R_vec(idx,:);

% mean values
res_mean = struct;
Ia_mean = zeros(T+1,1); 
Is_mean = Ia_mean;
Iu_mean = Ia_mean;
Io_mean = Ia_mean;
E_mean = Ia_mean;
dI_mean = Ia_mean;
H_mean = Ia_mean;
C_mean = Ia_mean;
V_mean = Ia_mean;
N_mean = Ia_mean;
R_mean = Ia_mean;
S_mean = Ia_mean;
for t = 1:T+1
    S_mean(t) = mean(S_vec(:,t));
    E_mean(t) = mean(E_vec(:,t));
    Ia_mean(t) = mean(Ia_vec(:,t));
    Is_mean(t) = mean(Is_vec(:,t));
    Io_mean(t) = mean(Io_vec(:,t));
    Iu_mean(t) = mean(Iu_vec(:,t));
    dI_mean(t) = mean(dI_in_vec(:,t));
    N_mean(t) = mean(N_vec(:,t));
    V_mean(t) = mean(V_vec(:,t));
    H_mean(t) = mean(H_vec(:,t));
    C_mean(t) = mean(C_vec(:,t));
    R_mean(t) = mean(R_vec(:,t));
end
% outputs
res_mean.Inf_asymp = Ia_mean;
res_mean.Inf_symp = Is_mean;
res_mean.Inf_obs = Io_mean;
res_mean.Inf_unobs = Iu_mean;
res_mean.dInf = dI_mean;
res_mean.Hosp_normal = N_mean;
res_mean.Hosp_icu = C_mean;
res_mean.Hosp_vent = V_mean;
res_mean.Hosp = H_mean;
res_mean.Rec = R_mean;
res_mean.Sus = S_mean;
res_mean.Exp = E_mean;

% quantiles & distribution
M = length(q_vec);
res_q.Inf_asymp = get_quant(Ia_vec);
res_q.Inf_symp = get_quant(Is_vec);
res_q.Inf_obs = get_quant(Io_vec);
res_q.Inf_unobs = get_quant(Iu_vec);
res_q.dInf = get_quant(dI_in_vec);
res_q.Hosp_normal = get_quant(N_vec);
res_q.Hosp_icu = get_quant(C_vec);
res_q.Hosp_vent = get_quant(V_vec);
res_q.Hosp = get_quant(H_vec);
res_q.Rec = get_quant(R_vec);
res_q.Sus = get_quant(S_vec);
res_q.Exp = get_quant(E_vec);

% helpers
    function [x] = get_rv(y)
        shape0 = y.mean*(y.std)^2; scale0 = 1/(y.std)^2;
        shape0_vec = shape0*ones(1*N,1);
        scale0_vec = scale0*ones(1*N,1);
        x = reshape(gamrnd(shape0_vec,scale0_vec),N,1);
    end

    function [x] = get_quant(vec)
        x = zeros(M,T);
        for j = 1:M
            x(j,:) = quantile(vec,q_vec(j),1);
        end
    end

    function [r] = get_kappa_res()
        dr = s.kappa_res_delta_0+restrictions.delta*min([0:T],restrictions.at)/restrictions.at; %#ok<NBRAK>
        r = (1+dr.^s.kappa_res_alpha).*s.beta_res;
    end

    function [m] = get_kappa_mob()
        m0 = interp1(mobility.x_grid,mobility.y_grid,mobility.values);
        ms = interp1(mobility.x_grid,mobility.y_grid,mobility.scale);
        m = m0/ms;
    end

end