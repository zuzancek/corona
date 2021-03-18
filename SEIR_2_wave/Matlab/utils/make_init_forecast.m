function []=make_init_forecast(s,data,dateFrom,dateTo)

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
% init: S, E, O, U
% R - reproduction number (effective)
% sigma_o, sigma_y

%% initialization
T = dateTo-dateFrom+1;
s.sim_num = s.pop_size;
N_o = ceil(s.sim_num.*s.dep_ratio_65);
N_y = s.sim_num-N_o;
alpha_o = s.alpha_i_o;
alpha_y = s.alpha_i_y;
alpha_s = s.alpha_s_o;
mu = s.mu;

% method_data = s.smoothing_method_data;

% transition times
T_lat_y = s.T_lat.mean;
T_inf_y = s.T_inf.mean; gamma_y = 1/T_inf_y;
T_lat_o = s.T_lat.mean;
T_inf_o = s.T_inf.mean; gamma_o = 1/T_inf_o;
T_pre_test_o = s.T_pre_test.mean;
T_pre_test_y = s.T_pre_test.mean;
T_inf_obs_o = s.T_inf_obs.mean;
T_inf_obs_y = s.T_inf_obs.mean;

% init values
S_o_ini = data.S_o;       S_y_ini = data.S_y;
E_o_ini = data.E_o;       E_y_ini = data.E_y;
O_o_ini = data.O_o;       O_y_ini = data.O_y;
U_o_ini = data.U_o;       U_y_ini = data.U_y;

Rt = data.Rt_avg; rt = double(Rt)+zeros(T,1); %
sigma_o = data.sigma_o_avg;
sigma_y = data.sigma_y_avg;

% arrays
U_o = zeros(T,1); O_o = U_o; E_o = O_o; S_o = O_o; X_o = O_o; V_o = O_o;
U_y = zeros(T,1); O_y = U_y; E_y = O_y; S_y = O_y; X_y = O_o; V_y = O_y;
% 
S_o(1) = S_o_ini;   S_y(1) = S_y_ini;
E_o(1) = E_o_ini;   E_y(1) = E_y_ini;
O_o(1) = O_o_ini;   O_y(1) = O_y_ini;
U_o(1) = U_o_ini;   U_y(1) = U_y_ini;

%% calculation
for t=2:T
    Z = (alpha_o*O_o(t-1)+alpha_s.*U_o(t-1))*gamma_o'/N_o+...
        (alpha_y*O_y(t-1)+U_y(t-1))*gamma_y'/N_y;
    S_o(t) = S_o(t-1)*(1-rt(t)*Z*mu);
    S_y(t) = S_y(t-1)*(1-rt(t)*Z);
    E_o(t) = E_o(t-1)*(1-1/T_lat_o)+S_o(t-1)*Z*rt(t)*mu;
    V_o(t) = E_o(t-1)/T_lat_o;
    E_y(t) = E_y(t-1)*(1-1/T_lat_y)+S_y(t-1)*Z*rt(t);
    V_y(t) = E_y(t-1)/T_lat_y;
    X_o(t) = sigma_o/T_pre_test_o*U_o(t-1);
    U_o(t) = U_o(t-1)-(1-sigma_o)*U_o(t-1)/T_inf_o-X_o(t)+V_o(t);
    X_y(t) = sigma_y/T_pre_test_y*U_y(t-1);
    U_y(t) = U_y(t-1)-(1-sigma_y)*U_y(t-1)/T_inf_y-X_y(t)+V_y(t);
    O_o(t) = O_o(t-1)*(1-1/T_inf_obs_o)+X_o(t);
    O_y(t) = O_y(t-1)*(1-1/T_inf_obs_y)+X_y(t);
end

%% data storage
p.S_o = tseries(dateFrom:dateTo,S_o);    p.S_y = tseries(dateFrom:dateTo,S_y);    p.S = p.S_o+p.S_y;
p.E_o = tseries(dateFrom:dateTo,E_o);    p.E_y = tseries(dateFrom:dateTo,E_y);    p.E = p.E_o+p.E_y;
p.O_o = tseries(dateFrom:dateTo,O_o);    p.O_y = tseries(dateFrom:dateTo,O_y);    p.O = p.O_o+p.O_y;
p.U_o = tseries(dateFrom:dateTo,U_o);    p.U_y = tseries(dateFrom:dateTo,U_y);    p.U = p.U_o+p.U_y;
p.X_o = tseries(dateFrom:dateTo,X_o);    p.X_y = tseries(dateFrom:dateTo,X_y);    p.X = p.X_o+p.X_y;
p.V_o = tseries(dateFrom:dateTo,V_o);    p.V_y = tseries(dateFrom:dateTo,V_y);    p.V = p.V_o+p.V_y;
p.I_o = p.O_o+p.U_o;                     p.I_y = p.O_y+p.U_y;                     p.I = p.I_o+p.I_y;

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

figure;
plot(p.V_o,'linewidth',2,'linestyle','--'); hold on;
plot(p.V_y,'linewidth',2,'linestyle','--');
plot(p.V,'linewidth',2,'color',0.5*[1 1 1],'linestyle',':');
plot(p.X_o,'linestyle','--','linewidth',2,'color','b');
plot(p.X_y,'linestyle','--','linewidth',2,'color','r');
plot(p.X,'linewidth',2,'color','k');
grid on;
title('New cases');
legend({'Unobserved - Old', 'Unobserved - Young', 'Unobserved - Total',...
    'Observed - Old', 'Observed - Young', 'Observed - Total'});

end