function [X,I,obs_ratio_adj] = adjust_infection_hospitals(x,h,s,dateFrom,dateTo,t0,t1)

%initialization
alpha_hr = s.alpha_hr;              % 7.73/100;
T_inf = s.T_inf.mean;               % 4.3;
T_hosp = s.T_hosp.mean;             % 4;
lambda = s.lambda;                  % 6.37/100
alpha_ih = lambda/T_hosp;
alpha_ir = 1/T_inf-alpha_ih;
T = dateTo-dateFrom+1;

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
    I(t) = d_I_H(t)/alpha_ih;
    d_I_R(t) = alpha_ir*I(t);
    I(t+1) = d_I_H(t+1)/alpha_ih;
    X(t) = I(t+1)-I(t)+d_I_H(t)+d_I_R(t);
end

% adjust series endpoints and get ratio
obs_ratio_adj = tseries(t0:t1,s.obs_ratio);
X(T-2:T) = X(T-3); X = [X(1);X(1);X]; % treat end-points
X0 = X;
X = tseries(dateFrom:dateTo+2,smooth_series(X,s.smooth_width,s.smooth_type,s.smooth_ends));
dI_data_real = resize(X,dateFrom:dateTo);
dI_data_reported = tseries(dateFrom:dateTo,dI_data);
delta = dI_data_reported./dI_data_real;

idx = find(dI_data_real<s.cases_min &dI_data_reported<s.cases_min & delta<1-s.ratio_threshold);
X0(idx) = dI_data_reported(idx);
dI_data_real(idx) = dI_data_reported(idx);
delta = dI_data_reported./dI_data_real;

obs_ratio_adj(dateFrom:dateTo) = smooth_series(delta*s.obs_ratio,s.smooth_width,s.smooth_type,s.smooth_ends);
XX = x.NewCases;
XX(dateFrom:dateTo) = tseries(dateFrom:dateTo,X0(1:T));
X = smooth_series(XX,s.smooth_width,s.smooth_type,s.smooth_ends);

end