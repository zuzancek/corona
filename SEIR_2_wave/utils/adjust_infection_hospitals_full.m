function [X,I,obs_ratio_adj,sa] = adjust_infection_hospitals_full(x,~,s,dateFrom,dateTo,t0,t1,sigma)


T = dateTo-dateFrom+1;

%
T_hosp = s.T_hosp.mean;
T_rec_norm = 7.9;
T_icu = 9.5-T_hosp;
T_rec_icu = 12.3-T_icu;
T_vent = 9-T_icu;
T_rec_vent = 15.5-(T_vent+T_icu);
T_death = 12-(T_vent+T_icu);
alpha_in = s.lambda./T_hosp.*ones(T,1);
alpha_ir = (1-lambda)/T_inf;

% alpha_vd = 7.62/100;
% alpha_cv = 9.86/100;
% alpha_nc = 7.61/100;
 
% % 
% % T_hosp = T_hosp0+zeros(T,1);
% % T_hosp(end) = T_hosp1; T_hosp(ceil(T/2)+1:end) = linspace(T_hosp0,T_hosp1,floor(T/2));
% % lambda = s.lambda;                  % 6.37/100

% definitions
D = smooth_series(x.Deaths(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
V = smooth_series(x.Ventilations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
C = smooth_series(x.ICU(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
H = smooth_series(x.Hospitalizations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
N = H-C-V;

alpha_vd = zeros(T,1); alpha_d = zeros(T,1); alpha_vr = zeros(T,1); 
alpha_cv = zeros(T,1); alpha_v = zeros(T,1); alpha_cr = zeros(T,1);
alpha_nc = zeros(T,1); alpha_c = zeros(T,1); alpha_nr = zeros(T,1);

d_V_R = zeros(T,1); 
d_C_V = zeros(T,1); d_C_R = zeros(T,1);
d_N_C = zeros(T,1); d_N_R = zeros(T,1);
d_I_N = zeros(T,1); d_I_R = zeros(T,1); I = zeros(T,1);
X = zeros(T,1);

% calculation
for t = 2:T
    alpha_vd(t) = (D(t)-D(t-1))/V(t-1);
    alpha_d(t) = T_death.*alpha_vd(t);
    alpha_vr(t) = (1-alpha_d(t))./T_rec_vent;
    d_V_R(t) = alpha_vr(t).*V(t-1);
    d_C_V(t) = V(t)-V(t-1)+d_V_D(t)+d_V_R(t);
    alpha_cv(t) = d_C_V(t)/C(t-1);
    alpha_v(t) = alpha_cv(t).*T_vent;
    alpha_cr(t) = (1-alpha_v(t))./T_rec_icu;
    d_C_R(t) = alpha_cr(t).*C(t-1);
    d_N_C(t) = C(t)-C(t-1)+d_C_V(t)+d_C_R(t);
    alpha_nc(t) = d_N_C(t)./N(t-1);
    alpha_c(t) = T_icu.*alpha_nc(t);
    alpha_nr(t) = (1-alpha_c(t))./T_rec_norm;
    d_N_R(t) = alpha_nr(t).*N(t-1);
    d_I_N(t) = N(t)-N(t-1)+d_N_C(t)+d_N_R(t);
    %
    I(t) = d_I_N(t)/alpha_in(t); 
    d_I_R(t) = I(t).*alpha_ir;
    X(t-1) = I(t)-I(t-1)+d_I_N(t-1)+d_I_R(t-1);
end

% adjust series endpoints and get ratio
obs_ratio_adj = tseries(t0:t1,s.obs_ratio);
dt = 2;
X(T-dt:T) = X(T-3); XX = X(T-3).*ones(dt,1); X = [XX;X]; X(1) = dI_data(1); % treat end-points
X0 = X; 
X = tseries(dateFrom:dateTo+dt,smooth_series(X,s.smooth_width,s.smooth_type,s.smooth_ends));
dI_data_real = resize(X,dateFrom:dateTo);
dI_data_reported = tseries(dateFrom:dateTo,dI_data);
delta = dI_data_reported./dI_data_real;

idx = find(dI_data_real<s.cases_min & dI_data_reported<s.cases_min & delta<1-s.ratio_threshold);
X0(idx-dateFrom+1) = dI_data_reported(idx); X0(1:min(idx)-dateFrom+1) = dI_data_reported(dateFrom:min(idx));
dI_data_real(idx) = dI_data_reported(idx);dI_data_real(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
delta = dI_data_reported./dI_data_real;

obs_ratio_adj(dateFrom:dateTo) = smooth_series(delta*s.obs_ratio,s.smooth_width,s.smooth_type,s.smooth_ends);
XX = x.NewCases;
XX(dateFrom:dateTo) = tseries(dateFrom:dateTo,X0(1:T));
X = smooth_series(XX,s.smooth_width,s.smooth_type,s.smooth_ends);

% [~,d] = max(resize(obs_ratio_adj,s.wave_2_from+5:dateTo));
symp_ratio_obs = 1-sigma(s.wave_2_from);
% [~,bp] = max((dI_data_reported(max(idx):dateTo)-dI_data_real(max(idx):dateTo)));
sa = struct;
sa.Xs = symp_ratio_obs.*X;
sa.Xa = X-sa.Xs;
sa.dIa_data_reported = dI_data_reported.*sigma(dateFrom:dateTo);
sa.dIs_data_reported = dI_data_reported-sa.dIa_data_reported;
sa.loss_a = sa.Xa-sa.dIa_data_reported;
sa.loss_s = sa.Xs-sa.dIs_data_reported;

end