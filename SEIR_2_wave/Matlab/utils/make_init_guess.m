function [p]=make_init_guess(s,data,dateFrom,dateTo)

%% scheme
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

%% initialization
T = dateTo-dateFrom+1;
s.sim_num = s.pop_size;
N_o = ceil(s.sim_num.*s.dep_ratio_65);
N_y = s.sim_num-N_o;
alpha_o = s.alpha_i_o;
alpha_y = s.alpha_i_y;
alpha_s = s.alpha_s_o;
mu = s.mu;

method_data = s.smoothing_method_data;

% transition times
T_lat_y = s.T_lat.mean;
T_inf_y = s.T_inf.mean; gamma_y = 1/T_inf_y;
T_lat_o = s.T_lat.mean;
T_inf_o = s.T_inf.mean; gamma_o = 1/T_inf_o;

% inputs
rho = remove_nan(data.rho,dateFrom,dateTo);
scale = s.sim_num/s.pop_size;
X_obs = scale*method_data(data.X_obs);
AC = scale*method_data(data.AC);
TC = scale*method_data(data.TC);
Rt = data.Rt; rt = double(Rt);

% new cases, infectious (observed/unoberved)
X_o_obs = X_obs.*rho;                           x_o_obs = double(X_o_obs);
X_y_obs = X_obs-X_o_obs;                        x_y_obs = double(X_y_obs);
X_y_unobs = X_y_obs./s.obs_ratio;               x_y_unobs = double(X_y_unobs);
X_o_unobs = s.alpha_s_o*X_o_obs./s.obs_ratio;   x_o_unobs = double(X_o_unobs);
x_o = x_o_obs+x_o_unobs;
x_y = x_y_obs+x_y_unobs;

U_o = zeros(T,1); O_o = U_o; E_o = O_o;
U_y = zeros(T,1); O_y = U_y; E_y = O_y;
sigma_o = s.obs_ratio*(1+0*double(rho)*(N_y+N_o)/N_o);
sigma_y = s.obs_ratio*(1+0*(1-double(rho))*(N_y+N_o)/N_y);
O_o(1) = AC(dateFrom)*rho(dateFrom);            O_y(1) = AC(dateFrom)-O_o(1);
U_o(1) = O_o(1).*s.alpha_s_o./s.obs_ratio;      U_y(1) = O_y(1)./s.obs_ratio;
E_o(1) = x_o_obs(2)*T_lat_o/s.obs_ratio;        E_y(1) = x_y_obs(2)*T_lat_o/s.obs_ratio;
S_y = zeros(T,1);           S_y(1) = N_y-TC(dateFrom)*(1-rho(dateFrom))./s.obs_ratio;         
S_o = zeros(T,1);           S_o(1) = N_o-TC(dateFrom)*(rho(dateFrom))./s.obs_ratio;           

%% calculation
for t=2:T
    O_o(t) = O_o(t-1).*(1-gamma_o)+x_o_obs(t);
    E_o(t-1) = x_o_obs(t)*T_lat_o/sigma_o(t);
    x_o_unobs(t) = (1-sigma_o(t))/T_lat_o*E_o(t-1);
    O_y(t) = O_y(t-1).*(1-gamma_y)+x_y_obs(t);
    E_y(t-1) = x_y_obs(t)*T_lat_y/sigma_y(t);
    x_y_unobs(t) = (1-sigma_y(t))/T_lat_y*E_y(t-1);
    U_o(t) = U_o(t-1).*(1-gamma_o)+x_o_unobs(t);
    U_y(t) = U_y(t-1).*(1-gamma_y)+x_y_unobs(t);
    x_y(t) = x_y_obs(t)+x_y_unobs(t);
    x_o(t) = x_o_obs(t)+x_o_unobs(t);
    Z = (alpha_o*O_o(t-1)+alpha_s.*U_o(t-1))*gamma_o'/N_o+...
        (alpha_y*O_y(t-1)+U_y(t-1))*gamma_y'/N_y;
    S_o(t) = S_o(t-1)*(1-mu*rt(t)*Z);
    S_y(t) = S_y(t-1)*(1-rt(t)*Z);
    E_o(t) = E_o(t-1)+Z*S_o(t-1)*rt(t)*mu-x_o(t);
    E_y(t) = E_y(t-1)+Z*S_y(t-1)*rt(t)-x_y(t);
end

%% data storage
p.S_o = tseries(dateFrom:dateTo,S_o);    p.S_y = tseries(dateFrom:dateTo,S_y);    p.S = p.S_o+p.S_y;
p.E_o = tseries(dateFrom:dateTo,E_o);    p.E_y = tseries(dateFrom:dateTo,E_y);    p.E = p.E_o+p.E_y;
p.O_o = tseries(dateFrom:dateTo,O_o);    p.O_y = tseries(dateFrom:dateTo,O_y);    p.O = p.O_o+p.O_y;
p.U_o = tseries(dateFrom:dateTo,U_o);    p.U_y = tseries(dateFrom:dateTo,U_y);    p.U = p.U_o+p.U_y;
p.I_o = p.O_o+p.U_o;                     p.I_y = p.O_y+p.U_y;                     p.I = p.I_o+p.I_y;
p.sigma_o = tseries(dateFrom:dateTo,sigma_o);
p.sigma_y = tseries(dateFrom:dateTo,sigma_y);
p.Rt = Rt;
p.X_o_obs = X_o_obs;       p.X_y_obs = X_y_obs;         p.X_obs = p.X_o_obs+p.X_y_obs;
p.X_o_unobs = X_o_unobs;   p.X_y_unobs = X_y_unobs;     p.X_unobs = p.X_o_unobs+p.X_y_unobs;
p.X_o = p.X_o_obs+p.X_o_unobs;      p.X_y = p.X_y_obs+p.X_y_unobs;  p.X = p.X_o+p.X_y;

%% plotting
figure;
subplot(2,1,1)
plot(p.S_o,'linewidth',1); hold on;
plot(p.S_y,'linewidth',1);
plot(p.S,'k','linewidth',2);
grid on;
title('Suspectible (S)');
legend({'Old','Young','Total'});

subplot(2,1,2)
plot(p.E_o,'linewidth',1); hold on;
plot(p.E_y,'linewidth',1);
plot(p.E,'k','linewidth',2);
grid on;
title('Exposed (E)');
legend({'Old','Young','Total'});

figure;
hh1=plot(p.U_o,'linewidth',2,'linestyle',':'); hold on;
hh2=plot(p.U_y,'linewidth',2,'linestyle',':');
plot(p.U,'linewidth',2,'color',0.5*[1 1 1],'linestyle',':');
plot(p.O_o,'linestyle','--','linewidth',2,'color',hh1.Color);
plot(p.O_y,'linestyle','--','linewidth',2,'color',hh2.Color);
plot(p.O,'linewidth',2,'color',0.5*[1 1 1],'linestyle','--');
plot(p.I_o,'b','linewidth',3);
plot(p.I_y,'r','linewidth',3);
plot(p.I,'k','linewidth',3);
grid on;
title('Infectious (I)');
legend({'Unobserved - Old', 'Unobserved - Young', 'Unobserved - Total',...
    'Observed - Old', 'Observed - Young', 'Observed - Total',...
    'Old', 'Young', 'Total'});

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
