function [Rt,q_mat,res,x_mat,Rt_last] = estimate_Rt_SEIR_aug(inputs,s,do_quant,do_weight)

% structure of inputs:
% I0: initial number of observed infectious
% z: daily data of inflow of newly observed infections
% obs_ratio: daily data of observed ratio of new infections (wrt all daily
%       data on new infections)
% asymp_ratio: daily data of observed ratio of asymptomatic new infections (wrt all daily
%       observed data on new infections)

% initialization
z = inputs.z;
T = length(z);
I0 = inputs.I0;
lambda = s.lambda;
try
    rho = inputs.obs_ratio;
    assert(length(rho)>=T);
catch err %#ok<*NASGU>
    rho = zeros(T,1)+s.obs_ratio;
end
try
    alpha = inputs.asymp_ratio;
    assert(length(alpha)>=T);
catch err
    alpha = zeros(T,1)+(1-s.symp_ratio_obs);
end
k = s.T_inf_asymp.mean/(s.T_inf_symp.mean);
theta = (1-alpha);
% rho = rho*k;

N = s.sim_num;
pop_size = s.pop_size;

T_lat = s.T_lat;
T_pre = s.T_pre;
T_inf_asymp = s.T_inf_asymp;
T_inf_symp = s.T_inf_symp;
T_inf_symp0 = s.T_inf_symp; T_inf_symp0.mean = T_inf_symp.mean-T_pre.mean;
T_inf_unobs = s.T_inf_unobs;
T_inf_hosp = s.T_inf_hosp;
T_inf_asymp_vec = get_rv(T_inf_asymp);
T_inf_symp_vec = get_rv(T_inf_symp);
T_inf_symp0_vec = get_rv(T_inf_symp0);
T_inf_unobs_vec = get_rv(T_inf_unobs);
T_inf_hosp_vec = get_rv(T_inf_hosp);
T_lat_vec = get_rv(T_lat);
T_pre_vec = get_rv(T_pre);
gamma_lat = 1./T_lat_vec;
gamma_pre = 1./T_pre_vec;
gamma_asymp = 1./T_inf_asymp_vec;
gamma_symp = 1./T_inf_symp_vec;
gamma_symp0 = 1./T_inf_symp0_vec;
gamma_unobs = 1./T_inf_unobs_vec;
gamma_hosp = 1./T_inf_hosp_vec;

% set initial values
S_vec = zeros(N,T);       S_vec(:,1) = pop_size-I0;
E_vec = zeros(N,T);       E_vec(:,1) = 0; % first few periods are irrelevant
Ia_vec = zeros(N,T);      Ia_vec(:,1) = rho(1)*(1-alpha(1))*I0;
Is_vec = zeros(N,T);      Is_vec(:,1) = rho(1)*alpha(1)*I0;
Iu_vec = zeros(N,T);      Iu_vec(:,1) = (1-rho(1))*I0;
Rt_vec = zeros(N,T);
Is_in = 0*Is_vec; Ia_in = 0*Is_vec; Iu_in = 0*Is_vec; Io_in = 0*Is_vec;

Rt = zeros(T,1); 
Iat = Rt; It = Rt; Iot = Rt; Iut = Rt; Ist = Rt; Et = Rt; St = Rt; Xt = Rt; 
Isint = Rt; Iaint = Rt; Ioint = Rt; Iuint = Rt;
idx = ones(N,1);

% model
% S(t+1) = S(t)-F(t);
% E(t+1) = E(t)+F(t)-E(t)/T_lat;
% Iu(t+1) = Iu(t)+E(t)/T_lat-tau/T_test*Iu(t)-(1-tau)/T_inf_u*Iu(t);
% Io(t+1) = Io(t)+tau/T_test*Iu(t)-[sigma*lambda/T_hosp+(1-sigma*lambda)/T_inf_o]*Io(t);
%
% input data
% z(t) = daily inflow of observed newly detected cases
% z(t) = tau/T_test*Iu(t)
% 
% reprod.number R(t)
% F(t) = R(t)*S(t)/pop_size*(Iu(t)/(varsigma*T_inf_u)+Ia(t)/T_inf_a+
%   Is(t)[lambda/T_hosp+(1-lambda)/T_inf_s]

for t = 1:T-1
    % E_vec(:,t) = alpha(t).*z(t)./rho(t);
    E_vec(:,t) = (z(t)-theta(t)*Ia_vec(:,t).*gamma_pre).*T_lat_vec./rho(t);
    % theta(:,t) = (1-alpha(t)).*z(t).*T_pre_vec./Ia_vec(:,t);
    % Is_vec(:,t+1) = Is_vec(:,t)+(1-alpha(t)).*z(t)...
    Is_in(:,t) = theta(t)*Ia_vec(:,t).*gamma_pre;
    Is_vec(:,t+1) = Is_vec(:,t)+Is_in(:,t)...
        -(lambda.*gamma_hosp+(1-lambda).*gamma_symp0).*Is_vec(:,t);
    Ia_in(:,t) = rho(t)*E_vec(:,t).*gamma_lat;
    Ia_vec(:,t+1) = Ia_vec(:,t)+rho(t)*E_vec(:,t).*gamma_lat...
        -(theta(t).*gamma_pre+(1-theta(t)).*gamma_asymp).*Ia_vec(:,t);
    Iu_vec(:,t+1) = Iu_vec(:,t)+(1-rho(t))*E_vec(:,t).*gamma_lat...
        -gamma_unobs.*Iu_vec(:,t);
    Iu_in(:,t) = (1-rho(t))*E_vec(:,t).*gamma_lat;
    Io_in(:,t) = Ia_in(:,t);
    E_vec(:,t+1) = (z(t+1)-theta(t)*Ia_vec(:,t+1).*gamma_pre).*T_lat_vec./rho(t);
    % Ia_vec(:,t+1) = Ia_vec(:,t)+alpha(t).*z(t)...
    %     -(theta(:,t).*gamma_pre+(1-theta(:,t)).*gamma_asymp).*Ia_vec(:,t);
    % Iu_vec(:,t+1) = Iu_vec(:,t)+(1-rho(t))*E_vec(:,t).*gamma_lat...
    %     -gamma_unobs.*Iu_vec(:,t);
    % E_vec(:,t+1) = (z(t+1)-theta(:,t).*Ia_vec(:,t+1).*gamma_pre).*T_lat_vec./rho(t);
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
Is_in = Is_in(idx,:);
Ia_in = Ia_in(idx,:);
Iu_in = Iu_in(idx,:);
Io_in = Io_in(idx,:);
Iu_vec = Iu_vec(idx,:);
Io_vec = Ia_vec+Is_vec;
I_vec = Iu_vec+Io_vec;
E_vec = E_vec(idx,:);
X_vec = I_vec+E_vec;

for t = 1:T
    Rt(t) = mean(Rt_vec(:,t));
    Iat(t) = mean(Ia_vec(:,t));
    Ist(t) = mean(Is_vec(:,t));
    Isint(t) = mean(Is_in(:,t));
    Iaint(t) = mean(Ia_in(:,t));
    Iuint(t) = mean(Iu_in(:,t));
    Ioint(t) = mean(Io_in(:,t));
    Iut(t) = mean(Iu_vec(:,t));
    Iot(t) = mean(Io_vec(:,t));
    It(t) = mean(I_vec(:,t));
    St(t) = mean(S_vec(:,t));
    Et(t) = mean(E_vec(:,t));
    Xt(t) = mean(X_vec(:,t));
end
Rt = Rt(1:end-1);

if do_quant
    q_vec = s.quant;
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

res.Iaint = Iaint;
res.Isint = Isint;
res.Ioint = Ioint;
res.Iuint = Iuint;
res.Iat = Iat;
res.Ist = Ist;
res.Iut = Iut;
res.Iot = Iot;
res.It = It;
res.St = St;
res.Et = Et;
res.Xt = Xt;

if do_weight
    weights = s.pweight;
    last_num = length(weights);
    nn = size(Rt_vec(:,T),1);
    weights_mat = repmat(weights,nn,1);
    Rt_last = Rt_vec(:,T-last_num+1:T);
    Rt_last = sum(Rt_last.*weights_mat,2);
else
    Rt_last = Rt_vec(:,T-1);
end

    function [x] = get_rv(y)
        shape0 = y.mean*(y.std)^2; scale0 = 1/(y.std)^2;
        shape0_vec = shape0*ones(1*N,1);
        scale0_vec = scale0*ones(1*N,1);
        x = reshape(gamrnd(shape0_vec,scale0_vec),N,1);
    end
end