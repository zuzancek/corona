function [res_mean,res_quant] = train_SEIR(t0,T,inputs,s)

%% handle inputs
mobility = inputs.mobility;
restrictions = inputs.restrictions;
q_vec = s.quant;
pop_size = s.pop_size;
init = inputs.init;
Rt = inputs.Rt(t0,t0+T);
T_test = inputs.T_test;

% init values (observed)
St = init.S(t0);
Iot = init.Io(t0);
if isfield(init,'Iu') 
    Iut = init.Iu(t0);
    Iut_next = init.Iu(t0+1);
else
    Iut = Iot*(1-1/s.obs_ratio);
    Iut_next = init.Io(t0+1)*(1-1/s.obs_ratio);    
end
if isfield(init,'E')
    Et = init.E(t0);
else    
    Et = (Iut_next-Iut.*(1-1./s.T_.mean)).*s.T_inc.mean;
end

% setup
N = length(Rt);
S_vec = zeros(N,T+1);       S_vec(:,1) = St;
E_vec = zeros(N,T+1);       E_vec(:,1) = Et;
Iu_vec = zeros(N,T+1);      Iu_vec(:,1) = Iu;
Io_vec = zeros(N,T+1);      Io_vec(:,1) = Io;
Ia_vec = zeros(N,T+1);      Ia_vec(:,1) = alpha(1)*Io; 
Is_vec = zeros(N,T+1);      Is_vec(:,1) = (1-alpha(1))*Io; 
R_eff = zeros(N,T+1);       R_eff(:,1) = Rt;
dE_in_vec = zeros(N,T+1);   dE_in_vec(:,1) = Rt;
st = zeros(T+1,1);
it = st; iot = st; iut = st; ist = st; iat = st;

%% initialize
T_inc = get_rv(s.T_inc);
T_inf_obs0 = get_rv(s.T_inf_obs0);
T_inf_unobs = get_rv(s.T_inf_unobs);
T_inf_obs0_vec = get_rv(T_inf_obs0);
gamma_inc = 1./T_inc;
gamma_unobs = 1./T_inf_unobs;
gamma_test = 1./T_test;

idx = ones(N,1);
kappa_res = get_kappa_res();
kappa_mob = get_kappa_mob();
varsigma_unobs = s.self_isolation_effect*ones(T,1);
varsigma_obs = s.case_isolation_effect*ones(T,1);


for t=1:T
    R_eff(:,t+1) = Rt.*kappa_mob(t).*kappa_res(t);
    dE_in_vec(:,t) = R_eff(:,t).*S_vec(:,t)/pop_size.*...
        (varsigma_obs(t).*Io_vec(:,t).*((1-lambda)*gamma_symp+lambda*gamma_hosp)...
         +varsigma_unobs(t).*Iu_vec(:,t).*gamma_unobs/varsigma);
    S_vec(:,t+1) = S_vec(:,t)-dE_in_vec(:,t);     
    E_vec(:,t+1) = E_vec(:,t).*(1-gamma_inc)+dE_in_vec(:,t);
    Iu_vec(:,t+1) = Iu_vec(:,t)+gamma_inc.*E_vec(:,t)-tau(t).*Iu_vec(:,t).*gamma_test...
        -(1-tau(t)).*Iu_vec(:,t).*gamma_unobs;
    Io_vec(:,t+1) = Io_vec(:,t)+tau(t).*Iu_vec(:,t).*gamma_test...
        -(lambda(t).*alpha(t).*gamma_hosp+(1-lambda(t).*alpha(t))).*Io_vec(:,t);
    Ia_vec(:,t) = alpha(t).*Io_vec(:,t); Is_vec(:,t) = Io_vec(:,t)-Ia_vec(:,t);
    idx = idx & Io_vec(:,t+1)>=0 & Iu_vec(:,t+1)>=0 & E_vec(:,t+1)>=0;
end
idx = find(idx>0);
S_vec = S_vec(idx,:);
% E_vec = E_vec(idx,:);
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
res_quant.Iat = get_quant(Iat_vec);
res_quant.Ist = get_quant(Ist_vec);
res_quant.Iut = get_quant(Iut_vec);
res_quant.Iot = get_quant(Iot_vec);
res_quant.It = get_quant(It_vec);
res_quant.St = get_quant(St_vec);


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
        dr = s.kappa_res_delta_0+restrictions.delta*min([0:T],restrictions.at)/restrictions.at; %#ok<NBRAK>
        r = (1+dr.^s.kappa_res_alpha).*s.beta_res;
    end

    function [m] = get_kappa_mob()
        m0 = interp1(mobility.x_grid,mobility.y_grid,mobility.values);
        ms = interp1(mobility.x_grid,mobility.y_grid,mobility.scale);
        m = m0/ms;
    end

end