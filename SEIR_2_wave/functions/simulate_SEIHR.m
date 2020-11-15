function [res_mean,res_quant] = simulate_SEIHR(T,Rt,init,alpha_vec,t0,s)

% init values (observed)
It = init.I;
St = init.S;
Ht = init.H;
Dt = init.D;
Vt = init.V;

% initialize
T_inf_asymp = s.T_inf_asymp;
T_inf_symp = s.T_inf_symp;
T_inc = s.T_inc;
q_vec = s.quant;
pop_size = s.pop_size;
obs_ratio = s.obs_ratio;
symp_ratio_obs = s.symp_ratio_obs;
eta = obs_ratio*symp_ratio_obs;
lambda = s.lambda;

% setup
N = length(Rt);
T_inc_vec = get_rv(T_inc);
T_inf_asymp_vec = get_rv(T_inf_asymp);
T_inf_symp_vec = get_rv(T_inf_symp);
T_inf_novent_vec = get_rv(T_inf_novent);

gamma_asymp = 1./T_inf_asymp_vec;
gamma_symp = 1./T_inf_symp_vec;
gamma_novent = 1./T_novent_rec_vec;
zeta = 1./T_inf_novent_vec;
delta_asymp = 1./(T_inf_asymp_vec+T_inc_vec-s.T_overlay);
delta_symp = 1./(T_inf_symp_vec+T_inc_vec-s.T_overlay);
delta_symp_novent = 1./(T_inf_novent_vec+T_inc_vec-s.T_overlay);

S_vec = zeros(N,T+1);       S_vec(:,1) = St(t0);
Ia_vec = zeros(N,T+1);      Ia_vec(:,1) = (1-eta)*It(t0);
Is_vec = zeros(N,T+1);      Is_vec(:,1) = eta*It(t0);
Iobs_vec = zeros(N,T+1);    Iobs_vec(:,1) = Is_vec(:,1)+(1-symp_ratio_obs)*obs_ratio/(1-eta)*Ia_vec(:,1);
R_vec = zeros(N,T+1);
H_vec = zeros(N,T+1);       H_vec(:,1) = Ht(t0);
V_vec = zeros(N,T+1);       V_vec(:,1) = Vt(t0);
N_vec = zeros(N,T+1);       N_vec(:,1) = Ht(t0)-Vt(t0);
D_vec = zeros(N,T+1);       D_vec(:,1) = Dt(t0);
dI_in_vec = zeros(N,T+1);   dI_in_vec(:,1) = Rt;
R_eff = zeros(N,T+1);       R_eff(:,1) = Rt;

idx = ones(N,1);

% calculation
for t=1:T
    % epidemiological part
    R_eff(:,t+1) = Rt.*kappa_mob(t).*kappa_res(t);
    dI_in_vec(:,t+1) = R_eff(:,t+1).*S_vec(:,t)/pop_size.*(Ia_vec(:,t).*delta_asymp+Is_vec(:,t).*((1-lambda)*delta_symp+lambda*delta_symp_novent));
    S_vec(:,t+1) = S_vec(:,t)-dI_in_vec(:,t+1);
    Ia_vec(:,t+1) = (1-gamma_asymp).*Ia_vec(:,t)+(1-eta)*dI_in_vec(:,t);
    Is_vec(:,t+1) = (1-(1-lambda)*gamma_symp-lambda*zeta).*Is_vec(:,t)+eta*dI_in_vec(:,t);
    Iobs_vec(:,t+1) = Is_vec(:,t+1)+(1-symp_ratio_obs)*obs_ratio/(1-eta)*Ia_vec(:,t+1);
    % clinical part
    N_vec(:,t+1) = lambda*zeta.*Is_vec(:,t)+N_vec(:,t).*(omega_novent.*(1-psi_novent)+(1-omega_novent).*(1-gamma_novent-xi));
    V_vec(:,t+1) = (1-omega_novent).*xi.*N_vec(:,t)+V_vec(:,t).*(omega_vent.*(1-psi_vent)+(1-omega_vent).*(1-gamma_vent));
    H_vec(:,t+1) = N_vec(:,t+1)+V_vec(:,t+1);
    % final stages
    R_vec(:,t+1) = R_vec(:,t)+gamma_asymp.*Ia_vec(:,t)+gamma_symp.*(1-lambda).*Is_vec(:,t)+(1-omega_novent).*gamma_novent.*N_vec(:,t)+(1-omega_vent).*gamma_vent.*V_vec(:,t);
    D_vec(:,t+1) = D_vec(:,t)+omega_novent.*psi_novent.*N_vec(:,t)+omega_vent.*psi_vent.*V_vec(:,t);
    idx = idx & Ia_vec(:,t+1)>0 & Is_vec(:,t+1)>0;
end

% adjust for valid indices
idx = find(idx>0);
S_vec = S_vec(idx,:);
Ia_vec = Ia_vec(idx,:);
Is_vec = Is_vec(idx,:);
Iobs_vec = Iobs_vec(idx,:);
dI_in_vec = dI_in_vec(idx,:);
N_vec = N_vec(idx,:);
V_vec = V_vec(idx,:);
H_vec = H_vec(idx,:);
R_vec = R_vec(idx,:);

% mean values
res_mean = struct;
Ia_mean = zeros(T+1,1); 
Is_mean = Ia_mean;
Iobs_mean = Ia_mean;
dI_mean = Ia_mean;
H_mean = Ia_mean;
V_mean = Ia_mean;
N_mean = Ia_mean;
R_mean = Ia_mean;
for t = 1:T+1
    Ia_mean(t) = mean(Ia_vec(:,t));
    Is_mean(t) = mean(Is_vec(:,t));
    Iobs_mean(t) = mean(Iobs_vec(:,t));
    dI_mean(t) = mean(dI_in_vec(:,t));
    N_mean(t) = mean(N_vec(:,t));
    V_mean(t) = mean(V_vec(:,t));
    H_mean(t) = mean(H_vec(:,t));
    R_mean(t) = mean(R_vec(:,t));
end

res_mean.Inf_asymp = Ia_mean;
res_mean.Inf_symp = Is_mean;
res_mean.Inf_obs = Iobs_mean;
res_mean.dInf = dI_mean;
res_mean.Hosp_novent = N_mean;
res_mean.Hosp_vent = V_mean;
res_mean.Hosp = H_mean;
res_mean.Rec = R_mean;

M = length(q_vec);
It_quant = zeros(M,T+1);
dIt_quant = zeros(M,T+1);
St_quant = zeros(M,T+1);
Ft_quant = zeros(M,T+1);
res_quant = struct;
for j = 1:M
    dIt_quant(j,:) = quantile(dI_in_vec,q_vec(j),1);
    It_quant(j,:) = quantile(I_vec,q_vec(j),1);
    St_quant(j,:) = quantile(S_vec,q_vec(j),1);
    Ft_quant(j,:) = quantile(F_vec,q_vec(j),1);
end
res_quant.dIt = dIt_quant;
res_quant.It = It_quant;
res_quant.St = St_quant;
res_quant.Ft = Ft_quant;

    function [x] = get_rv(y)
        shape0 = y.mean*(y.std)^2; scale0 = 1/(y.std)^2;
        shape0_vec = shape0*ones(1*N,1);
        scale0_vec = scale0*ones(1*N,1);
        x = reshape(gamrnd(shape0_vec,scale0_vec),N,1);
    end
end