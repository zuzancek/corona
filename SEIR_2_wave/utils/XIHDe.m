function [out] = XIHDe(x,p,s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
I0 = init.I;    C0 = init.C;    V0 = init.V;   
N0 = init.H-V0-C0; D0 = init.D;

method_params = s.smoothing_method_params;
% rho = method_params(init.rho); 
rho = init.rho(dateFrom:dateTo);
varsigma = method_params(init.varsigma); varsigma = varsigma(dateFrom:dateTo);

% parameters
alpha_iry = p.alpha_iry;    alpha_iny = p.alpha_iny; 
alpha_iro = p.alpha_iro;    alpha_ino = p.alpha_ino; 
alpha_ncy = p.alpha_ncy;    alpha_nry = p.alpha_nry; 
alpha_nco = p.alpha_nco;    alpha_nro = p.alpha_nro; 
alpha_cvy = p.alpha_cvy;    alpha_cry = p.alpha_cry; 
alpha_cvo = p.alpha_cvo;    alpha_cro = p.alpha_cro; 
alpha_vdy = p.alpha_vdy;    alpha_vry = p.alpha_vry; 
alpha_vdo = p.alpha_vdo;    alpha_vro = p.alpha_vro; 

lambda_iny = p.lambda_iny; lambda_ino = p.lambda_ino;
lambda_ncy = p.lambda_ncy; lambda_nco = p.lambda_nco;
lambda_cvy = p.lambda_cvy; lambda_cvo = p.lambda_cvo;
old_share = s.old_share;

%%
% ******* Equations
% I(t+1) = I(t)+X(t)-I_N(t)-I_R(t);     
%       I_N(t) = alpha*I(t);   I_R(t) = zeta_i*I(t);
% N(t+1) = N(t)+I_N(t)-N_C(t)-N_R(t);
%       N_C(t) = beta*N(t);    N_R(t) = zeta_n*N(t);
% C(t+1) = C(t)+N_C(t)-C_V(t)-C_R(t);
%       C_V(t) = gamma*C(t);    C_R(t) = zeta_c*C(t);
% V(t+1) = V(t)+C_V(t)-V_D(t)-V_R(t);
%       V_D(t) = delta*V(t);    C_R(t) = zeta_v*V(t);
% D(t+1) = D(t)+V_D(t);

%%

% initialization
method_data = s.smoothing_method_data;
X = method_data(x.NewCases(dateFrom:dateTo));
X_o = X.*rho;   X_y = X-X_o;
D0 = method_data(D0(dateFrom:dateTo));
D_o = zeros(T,1); D_o(1) = D0(1)*varsigma(1); 
D_y = zeros(T,1); D_y(1) = D0(1)-D_o(1); 
N0 = method_data(N0(dateFrom:dateTo));
N_o = zeros(T,1); N_o(1) = N0(1)*lambda_ino*rho(1)/(lambda_ino*old_share+lambda_iny*(1-old_share)); 
N_y = zeros(T,1); N_y(1) = N0(1)-N_o(1); 
C0 = method_data(C0(dateFrom:dateTo));
C_o = zeros(T,1); C_o(1) = C0(1)*lambda_ino*lambda_nco*rho(1)/(lambda_ino*lambda_nco*old_share+lambda_iny*lambda_ncy*(1-old_share)); 
C_y = zeros(T,1); C_y(1) = C0(1)-C_o(1); 
V0 = method_data(V0(dateFrom:dateTo));
V_o = zeros(T,1); V_o(1) = V0(1)*lambda_ino*lambda_nco*lambda_cvo*rho(1)/(lambda_ino*lambda_nco*lambda_cvo*lambda_cvo*old_share+lambda_iny*lambda_ncy*lambda_cvy*(1-old_share)); 
V_y = zeros(T,1); V_y(1) = V0(1)-V_o(1); 
I0 = method_data(I0(dateFrom:dateTo));
I_o = zeros(T,1); I_o(1) = I0(1)*s.old_share;
I_y = zeros(T,1); I_y(1) = I0(1)-I_o(1);

d_I_R_y = zeros(T,1); d_I_N_y = d_I_R_y;    d_I_R_o = d_I_R_y;      d_I_N_o = d_I_R_y; 
d_N_R_y = d_I_R_y;    d_N_C_y = d_I_R_y;    d_N_R_o = d_I_R_y;      d_N_C_o = d_I_R_y;
d_C_R_y = d_I_R_y;    d_C_V_y = d_I_R_y;    d_C_R_o = d_I_R_y;      d_C_V_o = d_I_R_y;
d_V_R_y = d_I_R_y;    d_V_D_y = d_I_R_y;    d_V_R_o = d_I_R_y;      d_V_D_o = d_I_R_y;

% calculation
for t=1:T-1
    d_I_R_y(t) = alpha_iry(t).*I_y(t); d_I_N_y(t) = alpha_iny.*I_y(t);
    I_y(t+1) = I_y(t)+X_y(t)-d_I_R_y(t)-d_I_N_y(t);
    d_I_R_o(t) = alpha_iro(t).*I_o(t); d_I_N_o(t) = alpha_ino.*I_o(t);
    I_o(t+1) = I_o(t)+X_o(t)-d_I_R_o(t)-d_I_N_o(t);
    d_N_R_y(t) = alpha_nry(t).*N_y(t); d_N_C_y(t) = alpha_ncy.*N_y(t);
    N_y(t+1) = N_y(t)+d_I_N_y(t)-d_N_R_y(t)-d_N_C_y(t);
    d_N_R_o(t) = alpha_nro(t).*N_o(t); d_N_C_o(t) = alpha_nco.*N_o(t);
    N_o(t+1) = N_o(t)+d_I_N_o(t)-d_N_R_o(t)-d_N_C_o(t);
    d_C_R_y(t) = alpha_cry.*C_y(t); d_C_V_y(t) = alpha_cvy.*C_y(t);
    C_y(t+1) = C_y(t)+d_N_C_y(t)-d_C_R_y(t)-d_C_V_y(t);
    d_C_R_o(t) = alpha_cro.*C_o(t); d_C_V_o(t) = alpha_cvo.*C_o(t);
    C_o(t+1) = C_o(t)+d_N_C_o(t)-d_C_R_o(t)-d_C_V_o(t);
    d_V_R_y(t) = alpha_vry.*V_y(t); d_V_D_y(t) = alpha_vdy.*V_y(t);
    V_y(t+1) = V_o(t)+d_C_V_y(t)-d_V_R_y(t)-d_V_D_y(t);
    d_V_R_o(t) = alpha_vro.*V_o(t); d_V_D_o(t) = alpha_vdo.*V_o(t);
    V_o(t+1) = V_o(t)+d_C_V_o(t)-d_V_R_o(t)-d_V_D_o(t);
    D_y(t+1) = D_y(t)+d_V_D_y(t);
    D_o(t+1) = D_o(t)+d_V_D_o(t);
end

% store results
out = struct;
out.X_y = tseries(dateFrom:dateTo,X_y);
out.X_o = tseries(dateFrom:dateTo,X_o);
out.I_y = tseries(dateFrom:dateTo,I_y);
out.I_o = tseries(dateFrom:dateTo,I_o);
out.N_y = tseries(dateFrom:dateTo,N_y);
out.N_o = tseries(dateFrom:dateTo,N_o);
out.C_y = tseries(dateFrom:dateTo,C_y);
out.C_o = tseries(dateFrom:dateTo,C_o);
out.V_y = tseries(dateFrom:dateTo,V_y);
out.V_o = tseries(dateFrom:dateTo,V_o);
out.D_y = tseries(dateFrom:dateTo,D_y);
out.D_o = tseries(dateFrom:dateTo,D_o);
out.H_y = tseries(dateFrom:dateTo,N_y+C_y+V_y);
out.H_o = tseries(dateFrom:dateTo,N_o+C_o+V_o);
out.X = out.X_o+out.X_y;
out.I = out.I_o+out.I_y;
out.N = out.N_o+out.N_y;
out.C = out.C_o+out.C_y;
out.V = out.V_o+out.V_y;
out.D = out.D_o+out.D_y;
out.H = out.H_o+out.H_y;

end