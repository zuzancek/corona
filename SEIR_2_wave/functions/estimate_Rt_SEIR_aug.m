function [Rt,q_mat,res,x_mat,Rt_last] = estimate_Rt_SEIR_aug(z,I0,s,q_vec,varargin)

% initialization
T = length(z);
N = s.sim_num;
pop_size = s.pop_size;

T_lat = s.T_lat;
T_pre = s.T_pre;
T_inf_asymp = s.T_inf_asymp;
T_inf_symp = s.T_inf_symp;
T_inf_unobs = s.T_inf_unobs;
T_inf_hosp = s.T_inf_hosp;

rho = s.obs_ratio;
sigma = s.symp_ratio_obs;
theta = s.p_a_s;
lambda = s.lambda;

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

% set initial values
S_vec = zeros(N,T);       S_vec(:,1) = pop_size-I0;
E_vec = zeros(N,T);       E_vec(:,1) = 0; % first few periods are irrelevant
Ia_vec = zeros(N,T);      Ia_vec(:,1) = rho*(1-sigma)*I0;
Is_vec = zeros(N,T);      Is_vec(:,1) = rho*sigma*I0;
Iu_vec = zeros(N,T);      Iu_vec(:,1) = (1-rho)*I0;
Rt_vec = zeros(N,T);

Rt = zeros(T,1); 
Iat = Rt; It = Rt; Iot = Rt; Iut = Rt; Ist = Rt; Et = Rt; St = Rt; Xt = Rt;
idx = ones(N,1);

% model
% S(t+1) = S(t)-F(t);
% E(t+1) = E(t)+F(t)-E(t)/T_lat;
% Iu(t+1) = Iu(t)+(1-rho)*E(t)/T_lat-Iu(t)/T_inf_u;
% Ia(t+1) = Ia(t)+(1-rho)*E(t)/T_lat-[theta/T_pre-(1-theta)/T_inf_a]*Ia(t);
% Ia(t+1) = Ia(t)+rho*E(t)/T_lat-[theta/T_pre-(1-theta)/T_inf_a]*Ia(t);
% Is(t+1) = Is(t)+theta*Ia(t)/T_pre-[lambda/T_hosp-(1-lambda)/T_inf_s]*Is(t);
%
% input data
% z(t) = rho*E(t)/T_lat+theta*Ia(t)/T_pre
% z(t+1) = rho*E(t+1)/T_lat+theta*Ia(t+1)/T_pre
% 
% reprod.number R(t)
% F(t) = R(t)*S(t)/pop_size*(Iu(t)/(varsigma*T_inf_u)+Ia(t)/T_inf_a+
%   Is(t)[lambda/T_hosp+(1-lambda)/T_inf_s]

for t = 1:T-1
    E_vec(:,t) = (z(t)-theta*Ia_vec(:,t).*gamma_pre).*T_lat_vec./rho;
    Is_vec(:,t+1) = Is_vec(:,t)+theta*Ia_vec(:,t).*gamma_pre...
        -(lambda.*gamma_hosp+(1-lambda).*gamma_symp).*Is_vec(:,t);
    Ia_vec(:,t+1) = Ia_vec(:,t)+rho*E_vec(:,t).*gamma_lat...
        -(theta.*gamma_pre+(1-theta).*gamma_asymp).*Ia_vec(:,t);
    Iu_vec(:,t+1) = Iu_vec(:,t)+(1-rho)*E_vec(:,t).*gamma_lat...
        -gamma_unobs.*Iu_vec(:,t);
    E_vec(:,t+1) = (z(t+1)-theta*Ia_vec(:,t+1).*gamma_pre).*T_lat_vec./rho;
    F = E_vec(:,t+1)-(1-gamma_lat).*E_vec(:,t);
    S_vec(:,t+1) = S_vec(:,t)-F;
    I = gamma_unobs.*Iu_vec(:,t)+gamma_asymp.*Ia_vec(:,t)...
       +(lambda.*gamma_hosp+(1-lambda).*gamma_symp).*Is_vec(:,t);
    Rt_vec(:,t) = pop_size./S_vec(:,t).*F./I;
    idx = idx & Is_vec(:,t)>0 & Ia_vec(:,t)>0 & Iu_vec(:,t)> 0 & E_vec(:,t)>0;
end
idx = find(idx>0);
Rt_vec = Rt_vec(idx,:);
S_vec = S_vec(idx,:);
Ia_vec = Ia_vec(idx,:);
Is_vec = Is_vec(idx,:);
Iu_vec = Iu_vec(idx,:);
Io_vec = Ia_vec+Is_vec;
I_vec = Iu_vec+Io_vec;
E_vec = E_vec(idx,:);
X_vec = I_vec+E_vec;

if ~isempty(varargin)
    weights = varargin{1};
    last_num = length(weights);
    nn = size(Rt_vec(:,T),1);
    weights_mat = repmat(weights,nn,1);
    Rt_last = Rt_vec(:,T-last_num+1:T);
    Rt_last = sum(Rt_last.*weights_mat,2);
else
    Rt_last = Rt_vec(:,T-1);
end

for t = 1:T
    Rt(t) = mean(Rt_vec(:,t));
    Iat(t) = mean(Ia_vec(:,t));
    Ist(t) = mean(Is_vec(:,t));
    Iut(t) = mean(Iu_vec(:,t));
    Iot(t) = mean(Io_vec(:,t));
    It(t) = mean(I_vec(:,t));
    St(t) = mean(S_vec(:,t));
    Et(t) = mean(E_vec(:,t));
    Xt(t) = mean(X_vec(:,t));
end
Rt = Rt(1:end-1);

if (nargin)>3
    M = length(q_vec);
    q_mat = zeros(M,T);
    for j = 1:M
        q_mat(j,:) = quantile(Rt_vec,q_vec(j),1);
    end
    x_mat = zeros(M,T);
    for j = 1:M
        x_mat(j,:) = quantile(X_vec(:,1:end),q_vec(j),1);
    end
else
    q_mat = [];
    x_mat = [];
end

res.Iat = Iat;
res.Ist = Ist;
res.Iut = Iut;
res.Iot = Iot;
res.It = It;
res.St = St;
res.Et = Et;
res.Xt = Xt;

    function [x] = get_rv(y)
        shape0 = y.mean*(y.std)^2; scale0 = 1/(y.std)^2;
        shape0_vec = shape0*ones(1*N,1);
        scale0_vec = scale0*ones(1*N,1);
        x = reshape(gamrnd(shape0_vec,scale0_vec),N,1);
    end
end