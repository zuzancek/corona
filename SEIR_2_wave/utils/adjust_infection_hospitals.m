function [X,I,obs_ratio_adj,sa] = adjust_infection_hospitals(x,h,s,dateFrom,dateTo,t0,t1,sigma)

% sigma = ratio of asymptomatical cases
%initialization
alpha_hr = s.alpha_hr;              % 7.73/100;
T_inf = s.T_inf.mean;               % 4.3;
T_hosp0 = s.T_hosp.mean;            % 4;
T_hosp1 = s.T_hosp.mean-1;          % 4;
T = dateTo-dateFrom+1;

T_hosp = T_hosp0+zeros(T,1);
T_hosp(end) = T_hosp1; T_hosp(ceil(T/2)+1:end) = linspace(T_hosp0,T_hosp1,floor(T/2));
lambda = s.lambda;                  % 6.37/100
alpha_ih = lambda./T_hosp;
alpha_ir = (1-lambda)/T_inf;

% definitions
D = x.Deaths(dateFrom:dateTo);
D = smooth_series(D,s.smooth_width,s.smooth_type,s.smooth_ends);
d_H_D = D(2:end)-D(1:end-1);
d_H_D = smooth_series(d_H_D,s.smooth_width,s.smooth_type,s.smooth_ends);
H = h.Hospitalizations(dateFrom:dateTo);
H = smooth_series(H,s.smooth_width,s.smooth_type,s.smooth_ends);
d_H_R = alpha_hr.*H;
d_I_H = zeros(T-1,1); d_I_R = zeros(T-1,1);
I = zeros(T-1,1);  X = zeros(T-2,1);
dI_data = smooth_series(x.NewCases(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);

% calculation
for t = 1:T-2
    d_I_H(t) = H(t+1)-H(t)+d_H_R(t)+d_H_D(t);
    d_I_H(t+1) = H(t+2)-H(t+1)+d_H_R(t+1)+d_H_D(t+1);
    I(t) = d_I_H(t)/alpha_ih(t);
    d_I_R(t) = alpha_ir*I(t);
    I(t+1) = d_I_H(t+1)/alpha_ih(t);
    X(t) = I(t+1)-I(t)+d_I_H(t)+d_I_R(t);
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