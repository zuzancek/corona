function [res_mean,res_quant] = train_SEIR(time_interval,inputs,s)

%% handle inputs
dateFrom = time_interval.dateFrom;
dateTo = time_interval.dateTo;
T = dateTo-dateFrom+1;
mobility = inputs.mobility;
restrictions = inputs.restrictions;
q_vec = s.quant;
pop_size = s.pop_size;
init = inputs.init;
Rt = inputs.Rt;
N = s.sim_num;
T_test0 = inputs.T_test;
alpha = inputs.asymp_ratio;
tau = inputs.obs_ratio;
lambda = s.lambda;
T_inf_unobs = get_rv_I_unobs();

% init values (observed)
St = init.S(1);
Iot = init.Io(1);
if isfield(init,'Iu') 
    Iut = init.Iu(1);
    Iut_next = init.Iu(2);
elseif isfield (init,'I')
    Iut = init.I(1)-Iot;
    Iut_next = init.I(2)-init.Io(2);
else
    Iut = Iot*(1-1/inputs.obs_ratio(1));
    Iut_next = init.Io(2)*(1-1/inputs.obs_ratio(2));    
end
if isfield(init,'E')
    Et = init.E(1);
else    
    Et = (Iut_next-Iut.*(1-1./s.T_inf_unobs.mean(1))).*s.T_inc.mean;
end

% setup
S_vec = zeros(N,T+1);       S_vec(:,1) = St;
E_vec = zeros(N,T+1);       E_vec(:,1) = Et;
Iu_vec = zeros(N,T+1);      Iu_vec(:,1) = Iut;
Io_vec = zeros(N,T+1);      Io_vec(:,1) = Iot;
Ia_vec = zeros(N,T+1);      Ia_vec(:,1) = alpha(1)*Iot; 
Is_vec = zeros(N,T+1);      Is_vec(:,1) = (1-alpha(1))*Iot; 
R_eff = zeros(N,T+1);       R_eff(:,1) = Rt(:,1);
dE_in_vec = zeros(N,T+1);   
st = zeros(T+1,1);
it = st; iot = st; iut = st; ist = st; iat = st;

%% initialize
T_lat = get_rv(s.T_lat);
T_pre = get_rv(s.T_pre);
T_test = repmat(T_test0',N,1)+repmat(T_pre,1,T);
T_inf_obs0 = get_rv(s.T_inf_obs0);
T_inf_obs = T_inf_obs0+T_test;
T_hosp = get_rv(s.T_hosp);
gamma_lat = 1./T_lat;
gamma_unobs = 1./T_inf_unobs;
gamma_obs = 1./T_inf_obs;
gamma_obs0 = 1./T_inf_obs0;
gamma_test = 1./T_test;
gamma_hosp = 1./T_hosp;

idx = ones(N,1);
kappa_res = get_kappa_res();
kappa_mob = get_kappa_mob();
varsigma_unobs = s.self_isolation_effect*ones(T,1);
varsigma_obs = s.case_isolation_effect*ones(T,1);

for t=1:T
    R_eff(:,t) = Rt(t).*kappa_mob(t).*kappa_res(t);
    dE_in_vec(:,t) = R_eff(:,t).*S_vec(:,t)/pop_size.*...
        (varsigma_obs(t).*Io_vec(:,t).*((1-lambda)*gamma_obs(:,t)+lambda*gamma_hosp)...
         +varsigma_unobs(t).*Iu_vec(:,t).*gamma_unobs(:,t));
    S_vec(:,t+1) = S_vec(:,t)-dE_in_vec(:,t);     
    E_vec(:,t+1) = E_vec(:,t).*(1-gamma_lat)+dE_in_vec(:,t);
    Iu_vec(:,t+1) = Iu_vec(:,t)+gamma_lat.*E_vec(:,t)-tau(t).*Iu_vec(:,t).*gamma_test(:,t)...
        -(1-tau(t)).*Iu_vec(:,t).*gamma_unobs(:,t);
    Io_vec(:,t+1) = Io_vec(:,t)+tau(t).*Iu_vec(:,t).*gamma_test(:,t)...
        -(lambda.*alpha(t).*gamma_hosp+(1-lambda.*alpha(t)).*gamma_obs0).*Io_vec(:,t);
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

    function [x] = get_rv_I_unobs()
        L = size(s.T_inf_obs.mean,2);
        % L = size(T_inf_asymp.mean,2);
        % sh0 = T_inf_asymp.mean.*repmat(share_asymp,N,1)+repmat((1-share_asymp),N,1)*T_inf_symp.mean;
        n0 = floor(N*s.share_reas);
        sh0 = (s.T_inf_obs.mean).*ones(N,1);
        x0 = gamrnd(sh0.*s.T_inf_obs.std.^2,1/s.T_inf_obs.std^2);
        x0 = repmat(x0(1:n0,:),1,T);
        s.T_inf_unobs.d1.mean = sh0; s.T_inf_unobs.d1.std = s.T_inf_asymp.std; 
        n1 = N-n0;
        if n1>0
            x1 = rand(n1,T); a0 = sh0; a1 = 10; sc = 2.3;
            x1 = (power_law(x1,a0(n0+1:end,:),a1,sc));
            s.T_inf_unobs.d2.a0 = a0; s.T_inf_unobs.d2.a1 = a1;s.T_inf_unobs.d2.k = sc;
            s.T_inf_unobs.s0 = s.share_reas;
            x = [x0;x1];
        else
            x = x0;
        end
        s.T_inf_unobs.mean = mean(x); s.T_inf_unobs.std = std(x);
    end

end