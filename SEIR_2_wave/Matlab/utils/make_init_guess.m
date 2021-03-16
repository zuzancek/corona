function []=make_init_guess(s,data)

% scheme
% Sy' = Sy-R/Tinf*Sy*I/N
% So' = So-R/Tinf*So*mu*I/N
% Ey' = Ey+R/Tinf*Sy*I/N-Ey/Tlat
% Eo' = Eo+R/Tinf*So*mu*I/N-Eo/Tlat
% Uy' = Uy+(1-zeta_y)*Ey/Tlat-Uy/Tinf;
% Uo' = Uo+(1-zeta_o)*Eo/Tlat-Uo/Tinf;
% Oy' = Oy+zeta_y*E_y/Tlat-Oy/Tinf;
% Oo' = Oo+zeta_o*E_o/Tlat-Oo/Tinf;

% data
% Xobs - inflow of observed new cases
% rho - share of 65+ in new cases
% R - reproduction number (effective)


N = s.pop_size;
N_y = s.pop_size_y;
N_o = s.pop_size_o;

method_data = s.smoothing_method_data;

T_lat_y = get_rv(s.T_lat,N_y);
T_inc_y = get_rv(s.T_inc,N_y);
T_lat_o = get_rv(s.T_lat,N_o);
T_inc_o = get_rv(s.T_inc,N_o);

X_obs = data.Xobs;
X_o_obs = X_obs.*data.rho;
X_y_obs = X_obs-X_o_obs;
X_y_unobs = X_y_obs./s.obs_ratio;
X_o_unobs = s.alpha_s_o*X_o_obs./s.obs_ratio;
X_o = method_data(X_o_unobs+X_o_obs);
X_y = method_data(X_y_unobs+X_y_obs);

E_o = X_o.*T_lat_o;
E_y = X_y.*T_lat_y;
Z_o = extend(E_o(2:end)-E_o(1:end-1),1)+X_o;
Z_y = extend(E_y(2:end)-E_y(1:end-1),1)+X_y;


    function [x] = get_rv(y,n)
        shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
        L = length(shape0);
        shape0_vec = repmat(shape0,n,1);
        scale0_vec = scale0*ones(n,L);
        x = gamrnd(shape0_vec,scale0_vec);
    end

    function [y] = extend(x,t0)
        [xlen,xwid] = size(x);
        z = x(1,:)+zeros(xlen+t0,xwid);
        z(t0+1:end,:) = x;
        y = method_data(z);
    end

end
