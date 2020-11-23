function [Rt,q_mat,res,x_mat,Rt_last] = estimate_Rt_SEIR(inputs,s,do_quant,do_weight)

% structure of inputs:
% I0: initial number of observed infectious
% z: daily data of inflow of newly observed infections
% obs_ratio: daily data of observed ratio of new infections (wrt all daily
%       data on new infections)
% asymp_ratio: daily data of observed ratio of asymptomatic new infections (wrt all daily
%       observed data on new infections)
% T_test: testing delay (after symptoms onset)

% initialization
z = inputs.z;
T = length(z);
I0 = inputs.I0;
lambda = s.lambda;
try
    tau = inputs.obs_ratio;
    assert(length(tau)>=T);
    assert(0==1);
catch err %#ok<*NASGU>
    tau = zeros(T,1)+s.obs_ratio;
end
try
    alpha = inputs.asymp_ratio;
    assert(length(alpha)>=T);
    assert(0==1);
catch err
    alpha = zeros(T,1)+(1-s.symp_ratio_obs);
end
try
    T_test = reshape(inputs.T_test,1,[]);
    assert(length(T_test)>=T);
catch err
    T_test = zeros(1,T)+s.T_test0;
end
try
    N = inputs.sim_num;
    assert(N~=s.sim_num);
catch err
    N = s.sim_num;
end
pop_size = s.pop_size;

T_lat = s.T_lat;
T_pre = s.T_pre;
T_inf_asymp = s.T_inf_asymp;
T_inf_symp = s.T_inf_symp;
T_inf_hosp = s.T_inf_hosp;
m0 = T_pre.mean+T_test;
T_inf_symp0 = s.T_inf_symp; T_inf_symp0.mean = T_inf_symp.mean-m0;
T_inf_asymp0 = s.T_inf_asymp; T_inf_asymp0.mean = T_inf_asymp.mean-m0;
T_inf_symp.mean= T_inf_symp0.mean+m0;
T_inf_asymp.mean = T_inf_asymp0.mean+m0;
T_inf_hosp0 = s.T_inf_hosp; T_inf_hosp0.mean = T_inf_hosp.mean-m0;
T_inf_hosp.mean = T_inf_hosp0.mean+m0;
if N>1
    u = zeros(N,1);    
    % symptoms onset
    T_pre_vec = T_pre.mean*ones(N,1); %get_rv(T_pre);    
    % test taken
    T_test_vec = repmat(T_test,N,1)+repmat(T_pre_vec,1,T);
    % infectious periods
    T_inf_asymp0_vec = get_rv(T_inf_asymp0);
    T_inf_symp0_vec = get_rv(T_inf_symp0);
    T_inf_symp_vec = get_rv(T_inf_symp); %T_inf_symp0_vec+T_test_vec;
    T_inf_asymp_vec = get_rv(T_inf_asymp);%T_inf_asymp0_vec+T_test_vec;
    share_reas = 1; %s.share_reas; 
    share_asymp = 1-s.symp_ratio_obs;
    T_inf_unobs_vec = get_rv_I_unobs(share_reas,share_asymp); 
    % hospital admission    
    T_inf_hosp0_vec = get_rv(T_inf_hosp0);
    T_inf_hosp_vec = get_rv(T_inf_hosp); 
    T_lat_vec = get_rv(T_lat);
else    
    T_pre_vec = T_pre.mean;
    T_test_vec = repmat(T_test,N,1)+repmat(T_pre_vec,1,T);
    T_inf_asymp0_vec = T_inf_asymp.mean-m0;
    T_inf_asymp_vec = T_inf_asymp0_vec+m0;
    T_inf_symp0_vec = T_inf_symp.mean-m0;
    T_inf_symp_vec = T_inf_symp0_vec+m0;
    share_reas = s.share_reas; 
    share_asymp = 1-s.symp_ratio_obs;
    T_inf_unobs_vec = T_inf_asymp_vec.*share_asymp+(1-share_asymp)*T_inf_symp_vec;
    T_inf_hosp_vec = T_inf_hosp.mean;
    T_inf_hosp0_vec = T_inf_hosp.mean-m0;
    T_lat_vec = T_lat.mean;
end
gamma_lat = 1./T_lat_vec;
gamma_pre = 1./T_pre_vec;
gamma_asymp = 1./T_inf_asymp_vec;
gamma_asymp0 = 1./T_inf_asymp0_vec;
gamma_symp = 1./T_inf_symp_vec;
gamma_symp0 = 1./T_inf_symp0_vec;
gamma_unobs = 1./T_inf_unobs_vec;
gamma_hosp = 1./T_inf_hosp_vec;
gamma_hosp0 = 1./T_inf_hosp0_vec;

T_pre_pd = makedist('Gamma','a',s.T_pre.mean*s.T_pre.std^2,'b',1/s.T_pre.std^2);
sigma = 0*s.symp_ratio_obs*(1-cdf(T_pre_pd,T_test+s.T_pre.mean));
varsigma = s.self_isolation_effect*ones(T,1);

% set initial values
S_vec = zeros(N,T);       S_vec(:,1) = pop_size-I0;
E_vec = zeros(N,T);       E_vec(:,1) = 0; % first few periods are irrelevant
Ia_vec = zeros(N,T);      Ia_vec(:,1) = tau(1)*(1-alpha(1))*I0;
Is_vec = zeros(N,T);      Is_vec(:,1) = tau(1)*alpha(1)*I0;
Iu_vec = zeros(N,T);      Iu_vec(:,1) = (1-tau(1))*I0;
Rt_vec = zeros(N,T);
Is_in = 0*Is_vec; Ia_in = 0*Is_vec; Iu_in = 0*Is_vec; Io_in = 0*Is_vec;

Rt = zeros(T-1,1); 
Iat = Rt; It = Rt; Iot = Rt; Iut = Rt; Ist = Rt; Et = Rt; St = Rt; Xt = Rt; 
Isint = Rt; Iaint = Rt; Ioint = Rt; Iuint = Rt;
idx = ones(N,1);

% model
% S(t+1) = S(t)-F(t);
% E(t+1) = E(t)+F(t)-E(t)/T_lat;
% Iu(t+1) = Iu(t)+E(t)/T_lat-tau(t)/T_test(t)*Iu(t)-(1-tau(t))/T_inf_u*Iu(t);
% Ia(t+1) = Ia(t)+alpha*tau(t)/T_test*Iu(t)-[sigma(t)/T_pre+(1-sigma(t))/T_inf_a0]*Ia(t);
% Is(t+1) = Is(t)+(1-alpha)*tau(t)/T_test*Iu(t)+sigma*Is(t)/T_pre-[lambda/T_hosp-(1-lambda)/T_inf_s0]*Is(t)
%
% input data
% z(t) = daily inflow of observed newly detected cases
% z(t) = tau(t)/T_test*Iu(t) = z_a(t)+z_s(t);
% z(t+1) = tau(t+1)/T_test*Iu(t+1)
% 
% reprod.number R(t)
% F(t) = R(t)*S(t)/pop_size*(Iu(t)/(varsigma*T_inf_u)+Ia(t)/T_inf_a+
%   Is(t)[lambda/T_hosp+(1-lambda)/T_inf_s]

for t = 1:T-2
    Is_in(:,t) = (1-alpha(t)).*z(t);
    Ia_in(:,t) = alpha(t).*z(t);
    Io_in(:,t) = z(t);
    Iu_vec(:,t+1) = z(t+1).*T_test_vec(:,t+1)./tau(t+1);    
    %Ia_vec(:,t+1) = Ia_vec(:,t)+Ia_in(:,t)-(sigma(t).*gamma_pre+(1-sigma(t)).*gamma_asymp0(:,t)).*Ia_vec(:,t);    
    %Is_vec(:,t+1) = Is_vec(:,t)+Is_in(:,t)+sigma(t).*gamma_pre-(lambda.*gamma_hosp0(:,t)+(1-lambda).*gamma_symp0(:,t)).*Is_vec(:,t);
    E_vec(:,t) = T_lat_vec.*(Iu_vec(:,t+1)-Iu_vec(:,t).*(1-tau(t)./T_test_vec(:,t)-(1-tau(t)).*gamma_unobs(:,t)));
    Is_in(:,t+1) = (1-alpha(t+1)).*z(t+1);
    Ia_in(:,t+1) = alpha(t+1).*z(t+1);
    Io_in(:,t+1) = z(t+1);
    Iu_vec(:,t+2) = z(t+2).*T_test_vec(:,t+2)./tau(t+2);    
    %Ia_vec(:,t+2) = Ia_vec(:,t+1)+Ia_in(:,t+1)-(sigma(t+1).*gamma_pre+(1-sigma(t+1)).*gamma_asymp0(:,t+1)).*Ia_vec(:,t+1);    
    %Is_vec(:,t+2) = Is_vec(:,t+1)+Is_in(:,t+1)+sigma(t+1).*gamma_pre-(lambda.*gamma_hosp0(:,t+1)+(1-lambda).*gamma_symp0(:,t+1)).*Is_vec(:,t+1);
    E_vec(:,t+1) = T_lat_vec.*(Iu_vec(:,t+2)-Iu_vec(:,t+1).*(1-tau(t+1)./T_test_vec(:,t+1)-(1-tau(t+1)).*gamma_unobs(:,t+1)));
    F = E_vec(:,t+1)-E_vec(:,t).*(1-gamma_lat);
    S_vec(:,t+1) = S_vec(:,t)-F;
    I = gamma_unobs(:,t).*Iu_vec(:,t)+0*gamma_asymp(:,t).*Ia_vec(:,t)/2.5...
       +0*((1-lambda).*gamma_symp(:,t)).*Is_vec(:,t)/2.5; 
    Rt_vec(:,t) = pop_size./S_vec(:,t).*F./I;
%     idx = idx & Is_vec(:,t+1)>0 & Ia_vec(:,t+1)>0 & Iu_vec(:,t+1)> 0 & E_vec(:,t)>0;
%     idx = idx & Is_vec(:,t+2)>0 & Ia_vec(:,t+2)>0 & Iu_vec(:,t+2)> 0 & E_vec(:,t+1)>0;
    idx = idx & S_vec(:,t+1)>0 & S_vec(:,t)>0 & F>0 & Iu_vec(:,t+1)> 0 & Iu_vec(:,t+2)> 0;
    disp(length(find(idx>0)));
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

for t = 1:T-1
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
        shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
        L = length(shape0);
        shape0_vec = repmat(shape0,N,1);
        scale0_vec = scale0*ones(N,L);
        x = gamrnd(shape0_vec,scale0_vec);
    end

    function [x] = get_rv_I_unobs(share_reas,share_asymp)
        L = size(T_inf_asymp.mean,2);
        sh0 = T_inf_asymp.mean.*repmat(share_asymp,N,1)+repmat((1-share_asymp),N,1)*T_inf_symp.mean;
        n0 = floor(N*share_reas);
        x0 = gamrnd(sh0*T_inf_asymp.std^2,1/T_inf_asymp.std^2);
        x0 = x0(1:n0,:);
        s.T_inf_unobs.d1.mean = sh0; s.T_inf_unobs.d1.std = s.T_inf_asymp.std; 
        n1 = N-n0;
        if n1>0
            x1 = rand(n1,1); a0 = sh0-2*s.T_inf_asymp.std; a1 = 10; sc = 2.3;
            x1 = repmat(power_law(x1,a0,a1,sc),1,T);
            s.T_inf_unobs.d2.a0 = a0; s.T_inf_unobs.d2.a1 = a1;s.T_inf_unobs.d2.k = sc;
            s.T_inf_unobs.s0 = share_reas;
            x = [x0;x1];
        else
            x = x0;
        end
    end

end