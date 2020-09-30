function [z,znorm] = simulate_SEIR(T_inf,T_inc,T,N,g1,g2,w,obs_ratio,obs)

shape = T_inf.mean*(T_inf.std)^2; scale = 1/(T_inf.std)^2;
gamma_mat = 1./gamrnd(shape,scale,T,N);
shape = T_inc.mean*(T_inc.std)^2; scale = 1/(T_inc.std)^2;
delta_mat = 1./gamrnd(shape,scale,T,N);

Rt_mat = get_R0(T,N,g1,g2,w);
% obs: already adjusted for N
x_init = get_first_infections(obs,obs_ratio,N);
% indicator: 0 = noninfectious, 1 = infectious

y = zeros(4,N);
y(1,:) = 1-x_init; y(3,:) = x_init;
z = zeros(4,T);

for t = 1:T
    Rt_vec = Rt_mat(t,:); gamma_vec = gamma_mat(t,:); delta_vec = delta_mat(t,:);
    betta_vec = Rt_vec.*gamma_vec;
    y0 = y;
    dS_vec_out = betta_vec.*y0(1,:).*y0(3,:)/N; dS_vec = -dS_vec_out;
    dE_vec_out = delta_vec.*y0(2,:); dE_vec = dS_vec_out-dE_vec_out;
    dI_vec_out = gamma_vec.*y0(3,:); dI_vec = dE_vec_out-dI_vec_out;
    dR_vec = dI_vec_out;
    y(1,:) = y0(1,:)+dS_vec;    z(1,t+1) = sum(y(1,:));
    y(2,:) = y0(2,:)+dE_vec;    z(2,t+1) = sum(y(2,:));
    y(3,:) = y0(3,:)+dI_vec;    z(3,t+1) = sum(y(3,:));
    y(4,:) = y0(4,:)+dR_vec;    z(4,t+1) = sum(y(4,:));    
end
znorm = z/N;

end

