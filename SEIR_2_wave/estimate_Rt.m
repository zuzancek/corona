function [Rt,q_mat] = estimate_Rt(dI_inflow,I0,pop_size,T_rem,N,q_vec)

T = length(dI_inflow);
shape = T_rem.mean*(T_rem.std)^2; scale = 1/(T_rem.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
Trec_mat = reshape(gamrnd(shape_vec,scale_vec),N,1);

% set initial values
S_vec = zeros(N,T); S_vec(:,1) = pop_size-I0;
I_vec = zeros(N,T); I_vec(:,1) = I0;
Rt_vec = zeros(N,T); 
Rt = zeros(T,1);
idx = ones(N,1);

% S(t+1) = S(t)-R(t)*gamma*S(t)*I(t)/pop_size;
% I(t+1) = I(t)+R(t)*gamma*S(t)*I(t)/pop_size-gamma*I(t);
% dI_inflow(t) = R(t)*gamma*S(t)*I(t)/pop_size;
for t = 1:T
   % gamma_vec = gamma_mat(:,t);
   Trec_vec = Trec_mat;%(:,t);
   % Rt_vec(:,t) = pop_size.*dI_inflow(t)./(gamma_vec.*S_vec(:,t).*I_vec(:,t));
   Rt_vec(:,t) = pop_size.*dI_inflow(t).*Trec_vec./(S_vec(:,t).*I_vec(:,t));
   S_vec(:,t+1) = S_vec(:,t)-dI_inflow(t);
   I_vec(:,t+1) = I_vec(:,t).*(1-1./Trec_vec)+dI_inflow(t);%+dI_inflow(t);
   idx = idx & I_vec(:,t+1)>0;
end
idx = find(idx>0);
Rt_vec = Rt_vec(idx,:); %#ok<FNDSB>

for t = 1:T
    Rt(t) = mean(Rt_vec(:,t));
end


if (nargin)>5
    M = length(q_vec);
    q_mat = zeros(M,T);
    for j = 1:M
        q_mat(j,:) = quantile(Rt_vec,q_vec(j),1);
    end
else
    q_mat = [];
end


end