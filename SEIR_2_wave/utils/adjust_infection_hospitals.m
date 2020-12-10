function [X,Is,obs_ratio_adj,sa] = adjust_infection_hospitals(x,h,s,dateFrom,dateTo,t0,t1,sigma,T_hosp_data)

T_inf = s.SI.mean;              
T_symp = s.T_hosp.mean; 
T_hosp0 = T_hosp_data.init+T_symp; 
T_hosp1 = T_hosp_data.final+T_symp; 
T = dateTo-dateFrom+1;
sigma = sigma(dateFrom:dateTo);

H = smooth_series(h.Hospitalizations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
T_hosp = T_hosp0+zeros(T,1);
bp = T_hosp_data.bp-dateFrom;
T_hosp(end) = T_hosp1; T_hosp(bp:end) = linspace(T_hosp0,T_hosp1,T-bp+1);
T_hosp = smooth_series(T_hosp,s.smooth_width,s.smooth_type,s.smooth_ends);
lambda = s.lambda;                  
alpha_ih = lambda./T_hosp;
alpha_ir = (1-lambda)/T_inf;

% definitions
D = smooth_series(x.Deaths(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
d_H_D = smooth_series(D(2:end)-D(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
d_H = smooth_series(H(2:end)-H(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
H = smooth_series(H,s.smooth_width,s.smooth_type,s.smooth_ends);
omega_vent = 26.21/100;         T_rec_vent = 14.5;  x_vent = smooth_series(h.Ventilation(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
omega_icu = 4.92/100;         T_rec_icu = 12.3;   x_icu = smooth_series(h.ICU(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
omega_norm = 2.71/100;      T_rec_norm = 7.9;   x_norm = H-x_icu-x_vent;
T_rec_hosp = (T_rec_vent.*x_vent+T_rec_icu.*x_icu+T_rec_norm.*x_norm)./x_norm;
omega_hosp = (omega_vent.*x_vent+omega_icu.*x_icu+omega_norm.*x_norm)./x_norm;
alpha_hr = (1-omega_hosp)./T_rec_hosp;
d_H_R = alpha_hr.*H;

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
Xs = adjust_tail(Xs,3);
X = adjust_tail(X,3);
Is = adjust_tail(Is,2);

figure;plot(dI_data);hold on;plot(X);plot(Xs,'k')

obs_ratio_adj = tseries(t0:t1,s.obs_ratio); 
X = tseries(dateFrom:dateTo,X);
dI_data_real = X;
dI_data_reported = tseries(dateFrom:dateTo,dI_data);
delta = dI_data_reported./dI_data_real;

obs_ratio_adj(dateFrom:dateTo) = smooth_series(delta*s.obs_ratio,s.smooth_width,s.smooth_type,s.smooth_ends);
dI_data_s = dI_data_reported.*(1-sigma);
dI_data_a = dI_data_reported-dI_data_s;

sa = struct;
Xs = tseries(dateFrom:dateTo,Xs);
sa.Xs = Xs;
sa.Xa = X-Xs;
sa.dIa_data_reported = dI_data_a;
sa.dIs_data_reported = dI_data_s;
sa.loss_a = sa.Xa-sa.dIa_data_reported;
sa.loss_s = sa.Xs-sa.dIs_data_reported;

Is = tseries(dateFrom:dateTo,Is);

    function [x] = adjust_tail(x,k)
        dx = x(T-k)-x(T-k-1);
        x(T-k+1) = x(T-k)+2/3*dx;
        x(T-k+2) = x(T-k+1)+1/3*dx;
        for j=3:k
            x(T-k+j) = x(T-k+j-1)+1/3*1/(j-1)*dx;
        end        
    end
end