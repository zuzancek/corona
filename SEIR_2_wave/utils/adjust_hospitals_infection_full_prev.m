function [out] = adjust_hospitals_infection_full_prev(x,p,s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
I0 = init.I;    C0 = init.C;    V0 = init.V;   
N0 = init.H-V0-C0; D0 = init.D;

T_hosp = p.T_hosp;
T_rec_norm = p.T_rec_norm;      % 6.9;
T_icu = p.T_icu;                % 8-T_hosp;
T_rec_icu = p.T_rec_icu;        % 12.3-T_icu;
T_vent = p.T_vent;              % 5-T_icu;
T_rec_vent = p.T_rec_vent;      % 14.5-(T_vent+T_icu);
T_death = p.T_death;            % 12-(T_vent+T_icu);
T_inf = s.SI.mean;

% definitions
alpha_d = mean(p.alpha_d)*ones(T-1,1);    alpha_vd = alpha_d./T_death;    alpha_vr = (1-alpha_d)./T_rec_vent;
alpha_v = mean(p.alpha_v)*ones(T-1,1);    alpha_cv = alpha_v./T_vent;     alpha_cr = (1-alpha_v)./T_rec_icu;
alpha_c = mean(p.alpha_c)*ones(T-1,1);    alpha_nc = alpha_c./T_icu;      alpha_nr = (1-alpha_c)./T_rec_norm;
lambda = p.lambda;      alpha_in = lambda./T_hosp;      alpha_ir = (1-lambda)./T_inf.*ones(T-1,1);
X = smooth_series(x.NewCases(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
I = zeros(T,1); I(1) = I0; d_I_R = zeros(T-1,1); d_I_N = zeros(T,1);
N = zeros(T,1); N(1) = N0; d_N_R = zeros(T-1,1); d_N_C = zeros(T,1);
C = zeros(T,1); C(1) = C0; d_C_R = zeros(T-1,1); d_C_V = zeros(T,1);
V = zeros(T,1); C(1) = C0; d_V_R = zeros(T-1,1); d_V_D = zeros(T,1);
D = zeros(T,1); D(1) = D0;

% calculation
for t=1:T-1
    d_I_R(t) = alpha_ir(t).*I(t); d_I_N(t) = alpha_in(t).*I(t);
    I(t+1) = I(t)+X(t)-d_I_R(t)-d_I_N(t);
    d_N_R(t) = alpha_nr(t).*N(t); d_N_C(t) = alpha_nc(t).*N(t);
    N(t+1) = N(t)+d_I_N(t)-d_N_R(t)-d_N_C(t);
    d_C_R(t) = alpha_cr(t).*C(t); d_C_V(t) = alpha_cv(t).*C(t);
    C(t+1) = C(t)+d_N_C(t)-d_C_R(t)-d_C_V(t);
    d_V_R(t) = alpha_vr(t).*V(t); d_V_D(t) = alpha_vd(t).*V(t);
    V(t+1) = V(t)+d_C_V(t)-d_V_R(t)-d_V_D(t);
    D(t+1) = D(t)+d_V_D(t);
end

out = struct;
out.X = tseries(dateFrom:dateTo,X);
out.I = tseries(dateFrom:dateTo,I);
out.N = tseries(dateFrom:dateTo,N);
out.C = tseries(dateFrom:dateTo,C);
out.V = tseries(dateFrom:dateTo,V);
out.D = tseries(dateFrom:dateTo,D);
out.H = tseries(dateFrom:dateTo,N+C+V);

end