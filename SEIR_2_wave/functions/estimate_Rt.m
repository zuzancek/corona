function [Rt,q_mat,It,Xt,x_mat,Rt_last,St] = estimate_Rt(dI_inflow,I0,pop_size,T_rem,N,q_vec,varargin)

T = length(dI_inflow);
shape = T_rem.mean*(T_rem.std)^2; scale = 1/(T_rem.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
Trec_mat = reshape(gamrnd(shape_vec,scale_vec),N,1);

T_inf.mean = 2.9; T_inf.std = T_rem.std;
shapeE = T_inf.mean*(T_inf.std)^2; scaleE = 1/(T_inf.std)^2;
shapeE_vec = shapeE*ones(1*N,1);
scaleE_vec = scaleE*ones(1*N,1);
Tinf_mat = reshape(gamrnd(shapeE_vec,scaleE_vec),N,1); %#ok<*NASGU>
T_inc.mean = 5.1; T_inc.std = T_rem.std;
shapeI = T_inc.mean*(T_inc.std)^2; scaleI = 1/(T_inc.std)^2;
shapeI_vec = shapeI*ones(1*N,1);
scaleI_vec = scaleI*ones(1*N,1);
Tinc_mat = reshape(gamrnd(shapeI_vec,scaleI_vec),N,1);

% set initial values
S_vec = zeros(N,T); S_vec(:,1) = pop_size-I0;
I_vec = zeros(N,T); I_vec(:,1) = I0;
E_vec = zeros(N,T);
Rt_vec = zeros(N,T); 
X_vec = zeros(N,T);
dE_inflow = zeros(N,T);
dE_outflow = zeros(N,T);
Rt = zeros(T,1); It = Rt; Xt = Rt;Et = Rt;St = Rt;
idx = ones(N,1);
d = T_rem.mean+14;

% S(t+1) = S(t)-R(t)*gamma*S(t)*I(t)/pop_size;
% I(t+1) = I(t)+R(t)*gamma*S(t)*I(t)/pop_size-gamma*I(t);
% dI_inflow(t) = R(t)*gamma*S(t)*I(t)/pop_size;
for t = 1:T
   % gamma_vec = gamma_mat(:,t);
   Trec_vec = Trec_mat;%(:,t);
   % Rt_vec(:,t) = pop_size.*dI_inflow(t)./(gamma_vec.*S_vec(:,t).*I_vec(:,t));
   Rt_vec(:,t) = pop_size.*dI_inflow(t).*Trec_vec./(S_vec(:,t).*I_vec(:,t));
   S_vec(:,t+1) = S_vec(:,t)-dI_inflow(t);
   dE_inflow(:,t) = Trec_mat./T_inf.mean.*dI_inflow(t); dE_outflow(:,t) = E_vec(:,t)./T_inc.mean;
   E_vec(:,t+1) = E_vec(:,t)+dE_inflow(:,t)-dE_outflow(:,t);
   I_vec(:,t+1) = I_vec(:,t).*(1-1./Trec_vec)+dI_inflow(t);%+dI_inflow(t);
   X_vec(:,t+1) = X_vec(:,t).*(1-1/d)+dE_inflow(:,t);
   idx = idx & I_vec(:,t+1)>0;
end
idx = find(idx>0);
Rt_vec = Rt_vec(idx,:);
S_vec = S_vec(idx,:);
I_vec = I_vec(idx,:);
X_vec = X_vec(idx,:);
E_vec = E_vec(idx,:);

if ~isempty(varargin)
    weights = varargin{1};
    last_num = length(weights);
    nn = size(Rt_vec(:,T),1);
    weights_mat = repmat(weights,nn,1);
    Rt_last = Rt_vec(:,T-last_num+1:T);
    Rt_last = sum(Rt_last.*weights_mat,2);
else
    Rt_last = Rt_vec(:,T);
end

for t = 1:T
    Rt(t) = mean(Rt_vec(:,t));
    It(t) = mean(I_vec(:,t));
    St(t) = mean(S_vec(:,t));
    Xt(t) = mean(X_vec(:,t));
    Et(t) = mean(E_vec(:,t));
end

if (nargin)>5
    M = length(q_vec);
    q_mat = zeros(M,T);
    for j = 1:M
        q_mat(j,:) = quantile(Rt_vec,q_vec(j),1);
    end
    x_mat = zeros(M,T);
    for j = 1:M
        x_mat(j,:) = quantile(X_vec(:,2:end),q_vec(j),1);
    end
else
    q_mat = [];
    x_mat = [];
end

end