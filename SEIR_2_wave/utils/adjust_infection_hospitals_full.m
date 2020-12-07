function [X,I,obs_ratio_adj,sa] = adjust_infection_hospitals_full(x,h,s,dateFrom,dateTo,t0,t1,sigma)


T = dateTo-dateFrom+1;

%
T_hosp = s.T_hosp.mean;
T_rec_norm = 7.9;
T_icu = 9.5-T_hosp;
T_rec_icu = 12.3-T_icu;
T_vent = 9-T_icu;
T_rec_vent = 15.5-(T_vent+T_icu);
T_death = 12-(T_vent+T_icu);
alpha_in = s.lambda./T_hosp.*ones(T-1,1);
T_inf = s.T_inf.mean;
lambda = 6.37/100;
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
V = smooth_series(h.Ventilation(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
C = smooth_series(h.ICU(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
H = smooth_series(h.Hospitalizations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
N = H-C-V;

X = zeros(T,1);

d_V_D = smooth_series(D(2:end)-D(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
alpha_vd = d_V_D./V(1:end-1);
alpha_d = T_death.*alpha_vd;
alpha_vr = (1-alpha_d)./T_rec_vent;
d_V_R = alpha_vr.*V(1:end-1);
d_V = V(2:end)-V(1:end-1);
d_C_V = d_V+d_V_D+d_V_R;
alpha_cv = d_C_V./C(1:end-1);
alpha_v = alpha_cv.*T_vent;
alpha_cr = (1-alpha_v)./T_rec_icu;
d_C_R = alpha_cr.*C(1:end-1);
d_C = C(2:end)-C(1:end-1);
d_N_C = d_C+d_C_V+d_C_R;
alpha_nc = d_N_C./N(1:end-1);
alpha_c = T_icu.*alpha_nc;
alpha_nr = (1-alpha_c)./T_rec_norm;
d_N_R = alpha_nr.*N(1:end-1);
d_N = N(2:end)-N(1:end-1);
d_I_N = d_N+d_N_C+d_N_R;
I = d_I_N./alpha_in;
d_I_R = I.*alpha_ir;
d_I = I(2:end)+I(1:end-1);
X = d_I+d_I_N+d_I_R;

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