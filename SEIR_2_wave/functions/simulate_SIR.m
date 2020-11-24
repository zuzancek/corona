function [res_mean,res_quant] = simulate_SIR(T,Rt,It,st,alpha_vec,t0,s)

% initialize
T_SI = s.T_SI;
q_vec = s.quant;
pop_size = s.pop_size;
w = s.alpha_weight;

% setup
N = length(Rt);
shape = T_SI.mean*(T_SI.std)^2; scale = 1/(T_SI.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
Trec_mat = reshape(gamrnd(shape_vec,scale_vec),N,1);
gamma_vec = 1./Trec_mat;

S_vec = zeros(N,T+1);   S_vec(:,1) = st(t0);
I_vec = zeros(N,T+1);   I_vec(:,1) = It(t0);
dI_in_vec = zeros(N,T+1);
alpha_vec = double(resize(alpha_vec,t0:t0+T));
idx = ones(N,1);
it = zeros(T+1,1); st = zeros(T+1,1); ft = zeros(T+1,1);
F_vec = zeros(N,T+1);F_vec(:,1) = Rt;

% calculation
for t=1:T
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

end