function [res_mean,res_quant] = train_SEIR(T,inputs,s)

%% handle inputs
Rt = inputs.Rt;
mobility = inputs.mobility;
restrictions = inputs.restrictions;
q_vec = s.quant;
pop_size = s.pop_size;
init = inputs.init;

% init values (observed)
Iot = init.Io;
Iut = init.Iu;
Et = init.E;
St = init.S;
Ht = init.H;
Dt = init.D;
Vt = init.V;

% setup
N = length(Rt);

%%
% initialize
shape = T_SI.mean*(T_SI.std)^2; scale = 1/(T_SI.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
Trec_mat = reshape(gamrnd(shape_vec,scale_vec),N,1);
gamma_vec = 1./Trec_mat;

S_vec = zeros(N,T+1);   S_vec(:,1) = st(t0);
I_vec = zeros(N,T+1);   I_vec(:,1) = It(t0);
dI_in_vec = zeros(N,T+1);
alpha_vec = double(resize(alpha_vec,t0:t0+T));
it = zeros(T+1,1); st = zeros(T+1,1); ft = zeros(T+1,1);
F_vec = zeros(N,T+1);F_vec(:,1) = Rt;
%%

Io_vec = zeros(N,T+1);      Io_vec(:,1) = Rt;
Iu_vec = zeros(N,T+1);      Iu_vec(:,1) = Rt;
dE_in_vec = zeros(N,T+1);   dE_in_vec(:,1) = Rt;
R_eff = zeros(N,T+1);       R_eff(:,1) = Rt;

idx = ones(N,1);
kappa_res = get_kappa_res();
kappa_mob = get_kappa_mob();
% calculation
for t=1:T
    R_eff(:,t+1) = Rt.*kappa_mob(t).*kappa_res(t);
    dE_in_vec(:,t) = R_eff(:,t).*S_vec(:,t)/pop_size.*...
        (Io_vec(:,t).*((1-lambda)*gamma_symp+lambda*gamma_hosp)+Iu_vec(:,t).*gamma_unobs/varsigma);
    
    F_vec(:,t+1) = Rt.*(1+w*(alpha_vec(t)-1)).*gamma_vec.*I_vec(:,t)/pop_size;
    dI_in_vec(:,t) = Rt.*(1+w*(alpha_vec(t)-1)).*gamma_vec.*S_vec(:,t).*I_vec(:,t)/pop_size;
    S_vec(:,t+1) = S_vec(:,t)-dI_in_vec(:,t);
    I_vec(:,t+1) = (1-gamma_vec).*I_vec(:,t)+dI_in_vec(:,t);
    idx = idx & I_vec(:,t+1)>0;
end
idx = find(idx>0);
S_vec = S_vec(idx,:);
F_vec = F_vec(idx,:);
I_vec = I_vec(idx,:);
dI_in_vec = dI_in_vec(idx,:);

% mean values
res_mean = struct;
for t = 1:T+1
    it(t) = mean(I_vec(:,t));
    dit(t) = mean(dI_in_vec(:,t)); %#ok<AGROW>
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