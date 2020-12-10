function [Xs,Is,obs_ratio_adj,sa] = adjust_infection_hospitals(x,h,s,dateFrom,dateTo,t0,t1,sigma)

% sigma = ratio of asymptomatical cases
%initialization
% alpha_hr = s.alpha_hr;              % 7.73/100;
T_inf = 6.5; % s.T_inf.mean;               % 4.3;
T_symp = 5; %  symptoms onset
T_hosp0 = T_symp; % s.T_hosp.mean;            % 4;
T_hosp1 = 1+T_symp; % s.T_hosp.mean;            % 4;
T = dateTo-dateFrom+1;
sigma = sigma(dateFrom:dateTo);

H = smooth_series(h.Hospitalizations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
omega_vent = 26.21/100;         T_rec_vent = 14.5;  x_vent = smooth_series(h.Ventilation(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
omega_icu = 4.92/100;         T_rec_icu = 12.3;   x_icu = smooth_series(h.ICU(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
omega_norm = 2.71/100;      T_rec_norm = 7.9;   x_norm = H-x_icu-x_vent;
T_rec_hosp = (T_rec_vent.*x_vent+T_rec_icu.*x_icu+T_rec_norm.*x_norm)./x_norm;
omega_hosp = (omega_vent.*x_vent+omega_icu.*x_icu+omega_norm.*x_norm)./x_norm;
alpha_hr = (1-omega_hosp)./T_rec_hosp;

% T_hosp = T_hosp0+zeros(T,1);
% T_hosp(end) = T_hosp1; T_hosp(ceil(T/2)+1:end) = linspace(T_hosp0,T_hosp1,floor(T/2));
T_hosp = T_hosp0+zeros(T,1);
T_hosp(end) = T_hosp1; T_hosp(52:end) = linspace(T_hosp0,T_hosp1,T-52+1);
T_hosp = smooth_series(T_hosp,s.smooth_width,s.smooth_type,s.smooth_ends);
lambda = 0.0910;%s.lambda;                  % 6.37/100
alpha_ih = lambda./T_hosp;
alpha_ir = (1-lambda)/T_inf;

% definitions
D = smooth_series(x.Deaths(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
d_H_D = smooth_series(diff(D),s.smooth_width,s.smooth_type,s.smooth_ends);
d_H_R = alpha_hr.*H;
d_H = smooth_series(diff(H),s.smooth_width,s.smooth_type,s.smooth_ends);
d_Is_H = zeros(T-1,1); d_Is_R = zeros(T-1,1);
Is = zeros(T-1,1);  Xs = zeros(T-2,1); X = Xs;
dI_data = smooth_series(x.NewCases(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);

% calculation
for t = 1:T-2
    d_Is_H(t) = d_H(t)+d_H_R(t)+d_H_D(t);
    d_Is_H(t+1) = d_H(t+1)+d_H_R(t+1)+d_H_D(t+1);
    Is(t) = d_Is_H(t)/alpha_ih(t);
    d_Is_R(t) = alpha_ir*Is(t);
    Is(t+1) = d_Is_H(t+1)/alpha_ih(t);
    Xs(t) = Is(t+1)-Is(t)+d_Is_H(t)+d_Is_R(t);
    X(t) = Xs(t)./(1-sigma(t));
end

% adjust series endpoints and get ratio
close all; figure;plot(dI_data);hold on;plot(X);plot(Xs,'k')
obs_ratio_adj = tseries(t0:t1,s.obs_ratio);
dt = 2;
X(T-dt:T) = X(T-3); XX = X(T-3).*ones(dt,1); X = [XX;X]; X(1) = dI_data(1); % treat end-points
X0 = X; 
X = tseries(dateFrom:dateTo+dt,smooth_series(X,s.smooth_width,s.smooth_type,s.smooth_ends));
dI_data_real = resize(Xs,dateFrom:dateTo);
dI_data_reported = tseries(dateFrom:dateTo,dI_data);
delta = dI_data_reported./dI_data_real;

idx = find(dI_data_real<s.cases_min & dI_data_reported<s.cases_min & delta<1-s.ratio_threshold);
X0(idx-dateFrom+1) = dI_data_reported(idx); X0(1:min(idx)-dateFrom+1) = dI_data_reported(dateFrom:min(idx));
dI_data_real(idx) = dI_data_reported(idx);dI_data_real(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
delta = dI_data_reported./dI_data_real;

obs_ratio_adj(dateFrom:dateTo) = smooth_series(delta*s.obs_ratio,s.smooth_width,s.smooth_type,s.smooth_ends);
XX = x.NewCases;
XX(dateFrom:dateTo) = tseries(dateFrom:dateTo,X0(1:T));
Xs = smooth_series(XX,s.smooth_width,s.smooth_type,s.smooth_ends);

% [~,d] = max(resize(obs_ratio_adj,s.wave_2_from+5:dateTo));
symp_ratio_obs = 1-sigma(s.wave_2_from);
% [~,bp] = max((dI_data_reported(max(idx):dateTo)-dI_data_real(max(idx):dateTo)));
sa = struct;
sa.Xs = symp_ratio_obs.*Xs;
sa.Xa = Xs-sa.Xs;
sa.dIa_data_reported = dI_data_reported.*sigma(dateFrom:dateTo);
sa.dIs_data_reported = dI_data_reported-sa.dIa_data_reported;
sa.loss_a = sa.Xa-sa.dIa_data_reported;
sa.loss_s = sa.Xs-sa.dIs_data_reported;

end