function [X,I] = adjust_infection_hospitals(x,h,s,dateFrom,dateTo)

alpha_hr = 7.73/100;
T_inf = 4.3;
T_hosp = 4;
alpha_ih = 6.37/100*(1/T_hosp);
alpha_ir = 1/T_inf-alpha_ih;
T = dateTo-dateFrom+1;

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

for t=1:T-2
    d_I_H(t) = H(t+1)-H(t)+d_H_R(t)+d_H_D(t);
    d_I_H(t+1) = H(t+2)-H(t+1)+d_H_R(t+1)+d_H_D(t+1);
    I(t) = d_I_H(t)/alpha_ih;
    d_I_R(t) = alpha_ir*I(t);
    I(t+1) = d_I_H(t+1)/alpha_ih;
    X(t) = I(t+1)-I(t)+d_I_H(t)+d_I_R(t);
end

dx = X(T-3)-X(T-4);
X(T-2:T) = X(T-3);
X = tseries(dateFrom+2:dateTo+2,smooth_series(X,s.smooth_width,s.smooth_type,s.smooth_ends));
dI_data = tseries(dateFrom:dateTo,dI_data);
figure;
plot(X);
hold on
plot(dI_data);
end