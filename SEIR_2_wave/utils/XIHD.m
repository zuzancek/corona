function [out] = XIHD(x,p,s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
I0 = init.I;     
H0 = init.H; 
D0 = init.D;

method_params = s.smoothing_method_params;
try
    rho = p.rho;
catch err
    rho = method_params(init.rho); rho = rho(dateFrom:dateTo);
end
varsigma = method_params(init.varsigma); varsigma = varsigma(dateFrom:dateTo);
     
alpha_ihy = s.eta_y/s.T_hosp_y;        
alpha_iho = s.eta_o/s.T_hosp_o;  
alpha_iry = (1-s.eta_y)./p.T_sick_y; 
alpha_iro = (1-s.eta_o)./p.T_sick_o;
alpha_hdy = p.alpha_hdy; %s.omega_y/s.T_death_y;
alpha_hdo = p.alpha_hdo; % s.omega_o/s.T_death_o;
alpha_hry = p.alpha_hry; % (1-s.omega_y)./p.T_rec_y; 
alpha_hro = p.alpha_hro; %(1-s.omega_o)./p.T_rec_o; 
eta_y = s.eta_y;                   
eta_o = s.eta_o;           

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);     
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
% D(t+1) = D(t)+H_D(t);

% initialization
method_data = s.smoothing_method_data;
X = method_data(x.NewCases(dateFrom:dateTo));
X_o = X.*rho;   X_y = X-X_o;
D0 = method_data(D0(dateFrom:dateTo));
D_o = zeros(T,1); D_o(1) = D0(1)*varsigma(1); 
D_y = zeros(T,1); D_y(1) = D0(1)-D_o(1); 
H0 = method_data(H0(dateFrom:dateTo));
H_o = zeros(T,1); H_o(1) = H0(1)*eta_o*rho(1)/(eta_o*s.old_share+eta_y*(1-s.old_share));
H_y = zeros(T,1); H_y(1) = H0(1)-H_o(1);
I0 = method_data(I0(dateFrom:dateTo));
I_o = zeros(T,1); I_o(1) = I0(1)*s.old_share;
I_y = zeros(T,1); I_y(1) = I0(1)-I_o(1);
d_I_H_o = zeros(T,1);  d_I_H_y = d_I_H_o; d_I_R_o = d_I_H_o; d_I_R_y = d_I_H_o;
d_H_D_o = d_I_H_o; d_H_D_y = d_I_H_o; d_H_R_o = d_I_H_o; d_H_R_y = d_I_H_o;

% calculation
for t=1:T-1
    d_I_H_o(t) = alpha_iho.*I_o(t);     d_I_R_o(t) = alpha_iro(t).*I_o(t);
    I_o(t+1) = I_o(t)+X_o(t)-d_I_H_o(t)-d_I_R_o(t);
    d_I_H_y(t) = alpha_ihy.*I_y(t);  d_I_R_y(t) = alpha_iry(t).*I_y(t);
    I_y(t+1) = I_y(t)+X_y(t)-d_I_H_y(t)-d_I_R_y(t);
    d_H_D_o(t) = alpha_hdo.*H_o(t);  d_H_R_o(t) = alpha_hro.*H_o(t);
    H_o(t+1) = H_o(t)+d_I_H_o(t)-d_H_D_o(t)-d_H_R_o(t);
    d_H_D_y(t) = alpha_hdy.*H_y(t);  d_H_R_y(t) = alpha_hry.*H_y(t);
    H_y(t+1) = H_y(t)+d_I_H_y(t)-d_H_D_y(t)-d_H_R_y(t);
    D_o(t+1) = D_o(t)+d_H_D_o(t);
    D_y(t+1) = D_y(t)+d_H_D_y(t);
end

% store results
out = struct;
out.X_o = tseries(dateFrom:dateTo,X_o);
out.X_y = tseries(dateFrom:dateTo,X_y);
out.I_o = tseries(dateFrom:dateTo,I_o);
out.I_y = tseries(dateFrom:dateTo,I_y);
out.H_o = tseries(dateFrom:dateTo,H_o);
out.H_y = tseries(dateFrom:dateTo,H_y);
out.D_o = tseries(dateFrom:dateTo,D_o);
out.D_y = tseries(dateFrom:dateTo,D_y);
out.X = out.X_o+out.X_y;
out.I = out.I_o+out.I_y;
out.H = out.H_o+out.H_y;
out.D = out.D_o+out.D_y;

end