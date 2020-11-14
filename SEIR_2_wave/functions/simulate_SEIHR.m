function [res_mean,res_quant] = simulate_SEIHR(T,Rt,It,st,alpha_vec,t0,s)

% initialize
T_inf_asymp = s.T_inf_asymp;
T_inf_symp = s.T_inf_symp;
T_inc = s.T_inc;
q_vec = s.quant;
pop_size = s.pop_size;
obs_ratio = s.obs_ratio;
symp_ratio_obs = s.symp_ratio_obs;
eta = obs_ratio*symp_ratio_obs;
nu = s.nu;

% setup
N = length(Rt);
shape = T_inc.mean*(T_inc.std)^2; scale = 1/(T_inc.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
T_inc_vec = reshape(gamrnd(shape_vec,scale_vec),N,1);
shape = T_inf_asymp.mean*(T_inf_asymp.std)^2; scale = 1/(T_inf_asymp.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
T_inf_asymp_vec = reshape(gamrnd(shape_vec,scale_vec),N,1);
gamma_asymp = 1./T_inf_asymp_vec;
delta_asymp = 1./(T_inf_asymp_vec+T_inc_vec-s.T_overlay);
shape = T_inf_symp.mean*(T_inf_symp.std)^2; scale = 1/(T_inf_symp.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
T_inf_symp_vec = reshape(gamrnd(shape_vec,scale_vec),N,1);
gamma_symp = 1./T_inf_symp_vec;
shape = T_inc.mean*(T_inc.std)^2; scale = 1/(T_inc.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
T_inc_vec = reshape(gamrnd(shape_vec,scale_vec),N,1);
delta_symp = 1./(T_inf_symp_vec+T_inc_vec-s.T_overlay);
shape = T_hosp_rec.mean*(T_hosp_rec.std)^2; scale = 1/(T_hosp_rec.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
T_hosp_rec_vec = reshape(gamrnd(shape_vec,scale_vec),N,1);
gamma_hosp = 1./T_hosp_rec_vec;

S_vec = zeros(N,T+1);   S_vec(:,1) = st(t0);
Ia_vec = zeros(N,T+1);   Ia_vec(:,1) = (1-eta)*It(t0);
Is_vec = zeros(N,T+1);   Is_vec(:,1) = eta*It(t0);
Iobs_vec = zeros(N,T+1);   Iobs_vec(:,1) = Is_vec(:,1)+(1-symp_ratio_obs)*obs_ratio/(1-eta)*Ia_vec(:,1);
Ra_vec = zeros(N,T+1);
Rs_vec = zeros(N,T+1);
idx = ones(N,1);
dI_in_vec = zeros(N,T+1);dI_in_vec(:,1) = Rt;
R_eff = zeros(N,T+1);R_eff(:,1) = Rt;

% calculation
for t=1:T
    R_eff(:,t+1) = Rt.*kappa_mob(t).*kappa_res(t);
    dI_in_vec(:,t+1) = R_eff(:,t+1).*S_vec(:,t)/pop_size.*(Ia_vec(:,t).*delta_asymp+Is_vec(:,t).*delta_symp);
    S_vec(:,t+1) = S_vec(:,t)-dI_in_vec(:,t+1);
    Ia_vec(:,t+1) = (1-gamma_asymp).*Ia_vec(:,t)+(1-eta)*dI_in_vec(:,t);
    Is_vec(:,t+1) = (1-gamma_symp).*Is_vec(:,t)+eta*dI_in_vec(:,t);
    Iobs_vec(:,t+1) = Is_vec(:,t+1)+(1-symp_ratio_obs)*obs_ratio/(1-eta)*Ia_vec(:,t+1);
    H_vec(:,t+1) = (1-gamma_hosp).*H_vec(:,t)+nu.*gamma_symp.*Is_vec(:,t);
    Ra_vec(:,t+1) = Ra_vec(:,t)+gamma_asymp.*Ia_vec(:,t);
    Rs_vec(:,t+1) = Rs_vec(:,t)+gamma_symp.*(1-nu).*Is_vec(:,t);
    idx = idx & Ia_vec(:,t+1)>0 & Is_vec(:,t+1)>0;
end
idx = find(idx>0);
S_vec = S_vec(idx,:);
Ia_vec = Ia_vec(idx,:);
Is_vec = Is_vec(idx,:);
Iobs_vec = Iobs_vec(idx,:);
dI_in_vec = dI_in_vec(idx,:);

% mean values
res_mean = struct;
for t = 1:T+1
    it(t) = mean(I_vec(:,t));
    dit(t) = mean(dI_in_vec(:,t));
    st(t) = mean(S_vec(:,t));
    ft(t) = mean(F_vec(:,t));
end
res_mean.Ft = ft;
res_mean.It = it;
res_mean.dIt = dit;
res_mean.St = st;

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

end