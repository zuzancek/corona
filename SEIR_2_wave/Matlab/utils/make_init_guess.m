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
N_o = ceil(s.sim_num.*s.dep_ratio_65);
N_y = s.sim_num-N_o;
alpha_o = s.alpha_i_o;
alpha_y = s.alpha_i_y;
mu = s.alpha_s_o;

method_data = s.smoothing_method_data;

% transition times
T_lat_y = get_rv(s.T_lat,N_y);
T_inf_y = get_rv(s.T_inf,N_y);  gamma_y = 1./T_inf_y';
T_lat_o = get_rv(s.T_lat,N_o);
T_inf_o = get_rv(s.T_inf,N_o);  gamma_o = 1./T_inf_o';

% inputs
rho = remove_nan(data.rho,dateFrom,dateTo);
scale = s.sim_num/s.pop_size;
X_obs = scale*method_data(data.X_obs);
AC = scale*method_data(data.AC);
TC = scale*method_data(data.TC);

% new cases, infectious (observed/unoberved)
X_o_obs = X_obs.*rho;                           x_o_obs = double(X_o_obs);
X_y_obs = X_obs-X_o_obs;                        x_y_obs = double(X_y_obs);
X_y_unobs = X_y_obs./s.obs_ratio;               x_y_unobs = double(X_y_unobs);
X_o_unobs = s.alpha_s_o*X_o_obs./s.obs_ratio;   x_o_unobs = double(X_o_unobs);
x_o = x_o_obs+x_o_unobs;
x_y = x_y_obs+x_y_unobs;

U_o = zeros(T+1,N_o); O_o = U_o;
U_y = zeros(T+1,N_y); O_y = U_y; 
O_o(1,:) = AC(dateFrom)*rho(dateFrom);              O_y(1,:) = AC(dateFrom)-O_o(1);
U_o(1,:) = O_o(1,:).*s.alpha_s_o./s.obs_ratio;      U_y(1,:) = O_y(1,:)./s.obs_ratio;
idx_y = ones(1,N_y);
idx_o = ones(1,N_o);
for t=1:T
    O_o(t+1,idx_o) = O_o(t,idx_o).*(1-gamma_o(idx_o))+x_o_obs(t);
    O_y(t+1,idx_y) = O_y(t,idx_y).*(1-gamma_y(idx_y))+x_y_obs(t);
    U_o(t+1,idx_o) = U_o(t,idx_o).*(1-gamma_o(idx_o))+x_o_unobs(t);
    U_y(t+1,idx_y) = U_y(t,idx_y).*(1-gamma_y(idx_y))+x_y_unobs(t);
    idx_o = idx_o & O_o(t+1,:)>=0 & U_o(t+1,:)>=0;
    idx_y = idx_y & O_y(t+1,:)>=0 & U_y(t+1,:)>=0;
end

% exposed and newly exposed
E_o = T_lat_o*(x_o_unobs+x_o_obs)'; d_E_o = E_o(2:end,:)-E_o(1:end-1,:);
E_y = T_lat_y*(x_y_unobs+x_y_obs)'; d_E_y = E_y(2:end,:)-E_y(1:end-1,:);
F_o = d_E_o+x_o'; 
F_y = d_E_y+x_y';

% contact rate with infectious
Z = dot(alpha_o*O_o(1:end-1,:)+mu.*U_o(1:end-1,:),gamma_o)/N_o+...
    dot(alpha_y*O_o(1:end-1,:)+U_y(1:end-1,:),gamma_y)/N_y;

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

    function [x] = remove_nan(x,t0,t1)
        if isnan(x(t0))
            i0 = find(isnan(x(t0:t1)));
            if length(i0)>1
                di0 = i0(2:end)-i0(1:end-1);
                ii0 = find(di0>1,1);
                if isempty(ii0)
                    i0 = i0(end);
                else
                    i0 = i0(ii0)-1;
                end
            else
                i0=2;
            end
            x(t0) = x(t0+i0); 
        end
        if isnan(x(t1))
            i1 = find(isnan(x(t0:t1)));
            if length(i1)>1
                di1 = i1(2:end)-i1(1:end-1);
                i1 = i1(find(di1>1,1)-1);
                if isempty(i1)
                    i1 = di1(end);
                end
            else
                i1=t1-1;
            end
            x(t1) = x(i1);
        end
        x = interp(x);   
    end

end
