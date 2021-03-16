function [p]=make_init_guess(s,data,dateFrom,dateTo)

% scheme
% Sy' = Sy-Fy,          Fy = R/Tinf*Sy*Z
% So' = So-Fo,          Fo = R/Tinf*So*mu*Z
% Ey' = Ey+Fy-Ey/Tlat
% Eo' = Eo+Fo-Eo/Tlat
% Uy' = Uy+(1-zeta_y)*Ey/Tlat-Uy/Tinf;
% Uo' = Uo+(1-zeta_o)*Eo/Tlat-Uo/Tinf;
% Oy' = Oy+zeta_y*E_y/Tlat-Oy/Tinf;
% Oo' = Oo+zeta_o*E_o/Tlat-Oo/Tinf;
% Z = [(alpha_o*Oo+mu*Uo)/No+(alpha_y*Oy+Uy)/Ny]
%
% data
% Xobs - inflow of observed new cases
% rho - share of 65+ in new cases
% R - reproduction number (effective)

% initialization
T = dateTo-dateFrom+1;
N_y = s.pop_size_y;
N_o = s.pop_size_o;
alpha_o = s.alpha_i_o;
alpha_y = s.alpha_i_y;
mu = s.alpha_s_o;

method_data = s.smoothing_method_data;

% transition times
T_lat_y = get_rv(s.T_lat,N_y);
T_inf_y = get_rv(s.T_inf,N_y);  gamma_y = 1./T_inf_y;
T_lat_o = get_rv(s.T_lat,N_o);
T_inf_o = get_rv(s.T_inf,N_o);  gamma_o = 1./T_inf_o;

% inputs
rho = method_data(data.rho);
X_obs = method_data(data.Xobs);
AC = method_data(data.AC);
TC = method_data(data.TC);

% new cases, infectious (observed/unoberved)
X_o_obs = X_obs.*rho;
X_y_obs = X_obs-X_o_obs;
X_y_unobs = X_y_obs./s.obs_ratio;
X_o_unobs = s.alpha_s_o*X_o_obs./s.obs_ratio;
U_o = zeros(T+1,N_o); O_o = U_o;
U_y = zeros(T+1,N_y); O_y = U_y; 
O_o(1,:) = AC(dateFrom)*rho(dateFrom);    O_y(1,:) = AC(dateFrom)-O_o(1);
U_o(1,:) = O_o.*s.alpha_s_o./s.obs_ratio; U_y(1,:) = O_y./s.obs_ratio;
for t=1:T
    O_o(t+1,:) = O_o(t,:).*(1-gamma_o)+X_o_obs(t);
    O_y(t+1,:) = O_y(t,:).*(1-gamma_y)+X_y_obs(t);
    U_o(t+1,:) = U_o(t,:).*(1-gamma_o)+X_o_unobs(t);
    U_y(t+1,:) = U_y(t,:).*(1-gamma_y)+X_y_unobs(t);
end

% exposed and newly exposed
E_o = (X_o_unobs+X_o_obs).*T_lat_o;
E_y = (X_y_unobs+X_y_obs).*T_lat_y;
F_o = extend(E_o(2:end)-E_o(1:end-1),1)+X_o; 
F_y = extend(E_y(2:end)-E_y(1:end-1),1)+X_y;

% contact rate with infectious
Z = (alpha_o*O_o(1:end-1,:)+mu.*U_o(1:end-1,:))/N_o+...
    (alpha_y*O_o(1:end-1,:)+U_y(1:end-1,:))/N_y;

S_y = zeros(T+1,N_y);           S_y(1,:) = N_y-TC(dateFrom)*(1-rho(dateFrom))./s.obs_ratio;         S_y_alt = S_y;
S_o = zeros(T+1,N_o);           S_o(1,:) = N_o-TC(dateFrom)*(rho(dateFrom))./s.obs_ratio;           S_o_alt = S_o;
d_S_o = 1-R.*gamma_o.*Z;        cd_S_o = cumprod(d_S_o,1);  S_o(2:end,:) = S_o(1,:).*cd_S_o;
d_S_y = 1-R.*gamma_y.*Z;        cd_S_y = cumprod(d_S_y,1);  S_y(2:end,:) = S_y(1,:).*cd_S_y;
cd_S_o_alt = -cumsum(F_o,1);    S_o_alt(2:end,:) = S_o(1,:)+cd_S_o_alt;
cd_S_y_alt = -cumsum(F_y,1);    S_y_alt(2:end,:) = S_y(1,:)+cd_S_y_alt;

p.S_o = S_o;    p.S_y = S_y;
p.E_o = E_o;    p.E_y = E_y;
p.O_o = O_o;    p.O_y = O_y;
p.U_o = U_o;    p.U_y = U_y;

figure;
plot(S_o,'linewidth',1); hold on;
plot(S_y,'linewidth',1);
plot(S_o_alt,'b--','linewidth',1);
plot(S_y_alt,'r--','linewidth',1);
grid on;
legend({'S_o','S_y','S_o alt','S_y alt'});

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
