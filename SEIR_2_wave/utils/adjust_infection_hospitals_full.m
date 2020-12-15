function [X,I,obs_ratio_adj,sa,p] = adjust_infection_hospitals_full(x,h,s,dateFrom,dateTo,t0,t1,sigma,omega,delay)

T = dateTo-dateFrom+1;

rho = s.old_share;

% delay in testing (gradual)
T_delay_0 = delay.v0;               T_delay_1 = delay.v1;               
T_delay_at = delay.at;
T_delay = zeros(T,1)+T_delay_0;     T_delay(T_delay_at-dateFrom:end) = T_delay_1;
T_delay = smooth_series(T_delay,s.smooth_width,s.smooth_type,s.smooth_ends);

T_inf = s.SI.mean;                  T_test = (2+1+s.T_pre.mean)+T_delay;
T_inf_y = T_inf-0.5;                T_inf_o = T_inf+1;
T_hosp_y = 7.02+T_test+1;           T_hosp_o = 3.24+T_test+1;         
lambda_y = 2.32/100;                lambda_o = 28.86/100;
alpha_h_y = lambda_y./T_hosp_y;     alpha_h_o = lambda_o./T_hosp_o;      alpha_h = rho*alpha_h_y+(1-rho)*alpha_h_o;
alpha_r_y = (1-lambda_y)/T_inf_y;   alpha_r_o = (1-lambda_o)/T_inf_o;      alpha_r = rho*alpha_r_o+(1-rho)*alpha_r_y;
T_death_y = 3.41;                   T_death_o = 4.59;
T_rec_y = 4.56;                     T_rec_o = 5.65;
omega_y = 5.15/100;                 omega_o = 37.16;
theta = rho/(1-rho)*lambda_o/lambda_y; 
theta = theta/(1+theta);
beta_d_y = omega_y/T_death_y;       beta_d_o = omega_o/T_death_o;        beta_d = theta*beta_d_o+(1-theta)*beta_d_y;

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);     
%       I_H(t) = alpha_h*I(t);   I_R(t) = alpha_r*I(t);
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
%       H_D(t) = beta_h*H(t);    H_R(t) = beta_r*H(t);
% D(t+1) = D(t)+H_D(t);

% initialization
dI_data = smooth_series(x.NewCases(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
D = smooth_series(x.Deaths(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
H = smooth_series(h.Hospitalizations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);

% calculation
H_D = smooth_series(D(2:end)-D(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
beta_d_t = H_D./H(1:end-1); k_d_t = beta_d_t./beta_d;
beta_r_t = theta*(1-omega_y*k_d_t)/T_rec_y+(1-theta)*(1-omega_o*k_d_t)/T_rec_o;
% H_R = beta_r_t.*H(1:end-1);
I_H = H(2:end)-H(1:end-1).*(1-beta_d_t-beta_r_t);
I = I_H./alpha_h(1:T-1);
% I_R = alpha_r.*I;
X = I(2:end)-I(1:end-1).*(1-alpha_r-alpha_h(1:T-2));

% % calculation
% d_H_D = smooth_series(D(2:end)-D(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
% H_y = (alpha_o.*H(1:end-1)-d_H_D)./(alpha_o-alpha_y);
% H_o = H-H_y;
% H_D_y = alpha_y.*H_y;     
% H_D_o = alpha_o.*H_o;
% H_R_y = (1-omega_y)./T_rec_y.*H_y;
% H_R_o = (1-omega_o)./T_rec_o.*H_o;
% d_H_y = smooth_series(H_y(2:end)-H_y(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
% d_H_o = smooth_series(H_o(2:end)-H_o(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
% I_H_y = d_H_y+H_R_y(1:end-1)+H_D_y(1:end-1);
% I_H_o = d_H_o+H_R_o(1:end-1)+H_D_o(1:end-1);
% I_y = I_H_y./delta_y;
% I_o = I_H_o./delta_o;
% I_R_y = gamma_y.*I_y;
% I_R_o = gamma_o.*I_y;
% d_I_y = I_y(2:end)-I_y(1:end-1);
% d_I_o = I_o(2:end)-I_o(1:end-1);
% X_y = d_I_y+I_R_y(1:end-1)+I_H_y(1:end-1);
% X_o = d_I_o+I_R_o(1:end-1)+I_H_o(1:end-1);
% X = X_y+X_o;

% alpha_xn_y = 2.32/100; alpha_xn_o = 31.86/100;
% alpha_xn = alpha_xn_o.*omega(dateFrom:dateTo)+alpha_xn_y.*(1-omega(dateFrom:dateTo));
% alpha_xn = alpha_xn(1:end-1);
% 
% % definitions
% D = smooth_series(x.Deaths(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
% V = smooth_series(h.Ventilation(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
% C = smooth_series(h.ICU(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
% H = smooth_series(h.Hospitalizations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
% N = H-C-V;

% H_D = smooth_series(D(2:end)-D(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
% alpha_vd = H_D./V(1:end-1);
% alpha_d = T_death.*alpha_vd;
% alpha_vr = (1-alpha_d)./T_rec_vent;
% d_V_R = alpha_vr.*V(1:end-1);
% d_V = V(2:end)-V(1:end-1);
% d_C_V = d_V+H_D+d_V_R;
% alpha_cv = d_C_V./C(1:end-1);
% alpha_v = alpha_cv.*T_vent;
% alpha_cr = (1-alpha_v)./T_rec_icu;
% d_C_R = alpha_cr.*C(1:end-1);
% d_C = C(2:end)-C(1:end-1);
% d_N_C = d_C+d_C_V+d_C_R;
% alpha_nc = d_N_C./N(1:end-1);
% alpha_c = T_icu.*alpha_nc;
% alpha_nr = (1-alpha_c)./T_rec_norm;
% d_N_R = alpha_nr.*N(1:end-1);
% d_N = N(2:end)-N(1:end-1);
% d_I_N = d_N+d_N_C+d_N_R;
% X = d_I_N./alpha_xn;

% I = d_I_N./alpha_in;
% d_I_R = I.*alpha_ir;
% d_I = I(2:end)-I(1:end-1);
% X = d_I+d_I_N(1:end-1)+d_I_R(1:end-1);

% adjust series endpoints and get ratio
obs_ratio_adj = tseries(t0:t1,s.obs_ratio);
X = adjust_tail(X,3);
X = tseries(dateFrom:dateTo,smooth_series(X,s.smooth_width,s.smooth_type,s.smooth_ends));
dI_data_real = resize(X,dateFrom:dateTo);
dI_data_reported = tseries(dateFrom:dateTo,dI_data);
delta = dI_data_reported./dI_data_real;

plot(X);hold on;% plot(dI_data_reported)
idx = find(dI_data_real<s.cases_min & dI_data_reported<s.cases_min & delta<1-s.ratio_threshold);
X(idx) = dI_data_reported(idx); X(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
dI_data_real(idx) = dI_data_reported(idx);dI_data_real(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
delta = dI_data_reported./dI_data_real;

obs_ratio_adj(dateFrom:dateTo) = smooth_series(delta*s.obs_ratio,s.smooth_width,s.smooth_type,s.smooth_ends);
XX = x.NewCases;
XX(dateFrom:dateTo) = X;
X = smooth_series(XX,s.smooth_width,s.smooth_type,s.smooth_ends);

sa = struct;
sa.Xs = (1-sigma(dateFrom)).*X;%
sa.Xa = X-sa.Xs;
sa.dIa_data_reported = dI_data_reported.*sigma(dateFrom:dateTo);
sa.dIs_data_reported = dI_data_reported-sa.dIa_data_reported;
sa.loss_a = sa.Xa-sa.dIa_data_reported;
sa.loss_s = sa.Xs-sa.dIs_data_reported; idx = find(sa.loss_s<0); sa.loss_s(idx) = 0; %#ok<FNDSB>

% store params & tseries
p.T_hosp = T_hosp;
p.T_rec_norm = T_rec_norm;
p.T_icu = T_icu; p.T_rec_icu = T_rec_icu;
p.T_vent = T_vent; p.T_rec_vent = T_rec_vent;
p.T_death = T_death; 
p.alpha_d = alpha_d; p.alpha_v = alpha_v; p.alpha_c = alpha_c;
p.lambda = lambda;

    function [x] = adjust_tail(x,k)
        dx = x(T-k)-x(T-k-1);
        x(T-k+1) = x(T-k)+2/3*dx;
        x(T-k+2) = x(T-k+1)+1/3*dx;
        for j=3:k
            x(T-k+j) = x(T-k+j-1)+1/3*1/(j-1)*dx;
        end        
    end

    function [x] = get_rv(y)
        shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
        L = length(shape0);
        shape0_vec = repmat(shape0,N,1);
        scale0_vec = scale0*ones(N,L);
        x = gamrnd(shape0_vec,scale0_vec);
    end

end