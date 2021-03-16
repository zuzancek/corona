function []=make_init_guess(s,data)

% scheme
% Sy' = Sy-R/Tinf*Sy*I/N
% So' = So-R/Tinf*So*mu*I/N
% Ey' = Ey+R/Tinf*Sy*I/N-Ey/Tlat
% Eo' = Eo+R/Tinf*So*mu*I/N-Eo/Tlat
% Iy' = Iy+Ey/Tlat-Iy/Tinf;
% Io' = Io+Eo/Tlat-Io/Tinf;

% data
% Xobs - inflow of observed new cases
% rho - share of 65+ in new cases
% R - reproduction number (effective)


N = s.pop_size;
N_y = s.pop_size_y;
N_o = s.pop_size_o;

T_lat_y = get_rv(s.T_lat,N_y);
T_inc_y = get_rv(s.T_inc,N_y);
T_lat_o = get_rv(s.T_lat,N_o);
T_inc_o = get_rv(s.T_inc,N_o);

X_obs = data.Xobs;
X_o_obs = X_obs.*data.rho;
X_y_obs = X_obs-X_o_obs;
X_y_unobs = X_y_obs./s.obs_ratio;
X_o_unobs = s.alpha_s_o*X_o_obs./s.obs_ratio;
X_o = X_o_unobs+X_o_obs;
X_y = X_y_unobs+X_y_obs;

E_o = X_o.*T_lat_o;
E_y = X_y.*T_lat_y;


function [x] = get_rv(y,n)
    shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
    L = length(shape0);
    shape0_vec = repmat(shape0,n,1);
    scale0_vec = scale0*ones(n,L);
    x = gamrnd(shape0_vec,scale0_vec);
end

end
