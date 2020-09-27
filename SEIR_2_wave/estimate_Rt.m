function [Rt,Rt_vec] = estimate_Rt(dI_inflow,S0,I0,pop_size,T_rem)

T = length(dI_inflow);
shape = T_rem.mean*(T_rem.std)^2; scale = 1/(T_rem.std)^2;
gamma_mat = 1./gamrnd(shape,scale,T,N);

% S_vec = [S0;dS_outflow_cum];
% set initial values
S_vec = zeros(T,N); S_vec(1,:) = S0;
I_vec = zeros(T,N); I_vec(1,:) = I0;
Rt_vec = zeros(T,N); 
Rt = zeros(T);

% S(t+1) = S(t)-R(t)*gamma*S(t)*I(t)/pop_size;
% I(t+1) = I(t)+R(t)*gamma*S(t)*I(t)/pop_size-gamma*I(t);
% dI_inflow(t) = R(t)*gamma*S(t)*I(t)/pop_size;
for t = 1:T
   gamma_vec = gamma_mat(t,:);
   Rt_vec(t,:) = pop_size.*dI_inflow(t)./(gamma_vec.*S_vec(t,:).*I_vec(t,:));
   S_vec(t+1,:) = S_vec(t,:)-dI_inflow(t);
   I_vec(t+1,:) = I_vec(t,:).*(1-gamma_vec)+dI_inflow(t);
   Rt(t) = mean(Rt_vec(t,:));
end