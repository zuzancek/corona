function [Rt,q_mat,It,Xt,x_mat,Rt_last,St,Et] = estimate_Rt_SEIR(dI_in,I0,pop_size,SI,N,q_vec,varargin)

T = length(dI_in);
T_inf.mean = 2.9; T_inf.std = SI.std;
shapeE = T_inf.mean*(T_inf.std)^2; scaleE = 1/(T_inf.std)^2;
shapeE_vec = shapeE*ones(1*N,1);
scaleE_vec = scaleE*ones(1*N,1);
T_inf_vec = reshape(gamrnd(shapeE_vec,scaleE_vec),N,1); %#ok<*NASGU>
T_inc.mean = 5.1; T_inc.std = SI.std;
shapeI = T_inc.mean*(T_inc.std)^2; scaleI = 1/(T_inc.std)^2;
shapeI_vec = shapeI*ones(1*N,1);
scaleI_vec = scaleI*ones(1*N,1);
T_inc_vec = reshape(gamrnd(shapeI_vec,scaleI_vec),N,1);

% set initial values
S_vec = zeros(N,T); S_vec(:,1) = pop_size-I0;
I_vec = zeros(N,T); I_vec(:,1) = I0;
E_vec = zeros(N,T);
Rt_vec = zeros(N,T);
dE_in = zeros(N,T);
Rt = zeros(T,1); It = Rt; Et = Rt; St = Rt; Xt = Rt;
idx = ones(N,1);

% S(t+1) = S(t)-R(t)/T_inf*S(t)*I(t)/pop_size;
% E(t+1) = E(t)+R(t)*gamma*S(t)*I(t)/pop_size-E(t)/T_inc;
% I(t+1) = I(t)+E(t)/T_inc-I(t)/T_inf;
% dI_in(t) = E(t)/T_inc;
% dE_in(t) = R(t)/T_inf*S(t)*I(t)/pop_size
for t = 1:T-1
    E_vec(:,t) = dI_in(t).*T_inc_vec;
    E_vec(:,t+1) = dI_in(t+1).*T_inc_vec;
    I_vec(:,t+1) = (1-1./T_inf_vec).*I_vec(:,t)+dI_in(t);
    dE_in(:,t) = E_vec(:,t+1)-(1-1./T_inc_vec).*E_vec(:,t);
    S_vec(:,t+1) = S_vec(:,t)-dE_in(:,t);
    Rt_vec(:,t) = pop_size.*dE_in(:,t).*T_inf_vec./(S_vec(:,t).*I_vec(:,t));
    idx = idx & I_vec(:,t+1)>0;
end
idx = find(idx>0);
Rt_vec = Rt_vec(idx,:);
S_vec = S_vec(idx,:);
I_vec = I_vec(idx,:);
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
    It(t) = mean(I_vec(:,t));
    St(t) = mean(S_vec(:,t));
    Et(t) = mean(E_vec(:,t));
    Xt(t) = mean(X_vec(:,t));
end
Rt = Rt(1:end-1);

if (nargin)>5
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

end