function [z,znorm] = simulate_SIR(T_rem,T,N,g1,g2,w,obs_ratio,obs)

shape = T_rem.mean*(T_rem.std)^2; scale = 1/(T_rem.std)^2;
gamma_mat = 1./gamrnd(shape,scale,T,N);

Rt_mat = get_R0(T,N,g1,g2,w);
% obs: already adjusted for N
x_init = get_first_infections(obs,obs_ratio,N);
% indicator: 0 = noninfectious, 1 = infectious

y = zeros(3,N);
y(1,:) = 1-x_init; y(2,:) = x_init;
z = zeros(3,T);

for t = 1:T
    Rt_vec = Rt_mat(t,:); gamma_vec = gamma_mat(t,:); 
    betta_vec = Rt_vec.*gamma_vec;
    y0 = y;
    dS_vec_out = betta_vec.*y0(1,:).*y0(2,:)/N; dS_vec = -dS_vec_out;
    dI_vec_out = gamma_vec.*y0(2,:); dI_vec = dS_vec_out-dI_vec_out;
    dR_vec = dI_vec_out;
    y(1,:) = y0(1,:)+dS_vec;    z(1,t+1) = sum(y(1,:));
    y(2,:) = y0(2,:)+dI_vec;    z(2,t+1) = sum(y(2,:));
    y(3,:) = y0(3,:)+dR_vec;    z(3,t+1) = sum(y(3,:));  
end
znorm = z/N;

end

