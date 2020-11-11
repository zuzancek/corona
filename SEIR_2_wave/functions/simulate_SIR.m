function [res_mean,res_quant] = simulate_SIR(T,Rt,It,St,alpha_vec,t0,s)

% initialize
T_rem = s.T_rem;
q_vec = s.quant;
pop_size = s.pop_size;
obs_ratio = s.obs_ratio;

% setup
N = length(Rt);
shape = T_rem.mean*(T_rem.std)^2; scale = 1/(T_rem.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
Trec_mat = reshape(gamrnd(shape_vec,scale_vec),N,1);
gamma_vec = 1./Trec_mat;

S_vec = zeros(N,T+1);   S_vec(:,1) = St(t0);
I_vec = zeros(N,T+1);   I_vec(:,1) = It(t0);
alpha_vec = alpha_vec(t0,end);
idx = ones(N,1);

% calculation
for t=1:T
    dI_in = Rt.*alpha_vec(t).*gamma.*S_vec(:,t).*I_vec(:,t)/pop_size;
    S_vec(:,t+1) = S_vec(:,t)-dI_in;
    I_vec(:,t+1) = (1-gamma_vec).*I_vec(:,t)+dI_in;
    idx = idx & I_vec(:,t+1)>0;
end
idx = find(idx>0);
S_vec = S_vec(idx,:);
I_vec = I_vec(idx,:);

% mean values
res_mean = struct;
for t = 1:T
    It(t) = mean(I_vec(:,t));
    St(t) = mean(S_vec(:,t));
end
res_mean.It = It;
res_mean.It_norm = It*obs_ratio;
res_mean.St = St;

M = length(q_vec);
It_quant = zeros(M,T);
St_quant = zeros(M,T);
res_quant = struct;
for j = 1:M
    It_quant(j,:) = quantile(I_vec,q_vec(j),1);
    St_quant(j,:) = quantile(S_vec,q_vec(j),1);
end
res_quant.It = It_quant;
res_quant.It_norm = It_quant*obs_ratio;
res_quant.St = St_quant;

end