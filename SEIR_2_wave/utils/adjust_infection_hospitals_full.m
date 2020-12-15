function [X,I,obs_ratio_adj,sa,p] = adjust_infection_hospitals_full(x,h,s,dateFrom,dateTo,t0,t1,sigma,omega,delay)

T = dateTo-dateFrom+1;

rho = s.old_share;

% delay in testing (gradual)
T_delay_0 = delay.v0;               T_delay_1 = delay.v1;               
T_delay_at = delay.at;
T_delay = zeros(T,1)+T_delay_0;     T_delay(T_delay_at-dateFrom:end) = T_delay_1;
T_delay = smooth_series(T_delay,s.smooth_width,s.smooth_type,s.smooth_ends);

T_inf = s.T_inf.mean;                   T_test = (1+s.T_pre.mean)+T_delay;
T_inf_y = T_inf;                        T_inf_o = T_inf+1;
T_hosp_y = 7.02+T_test;                 T_hosp_o = 3.24+T_test; 
dy = 0.5;                               do = 0.5;
T_rec_norm_y = 6.9-dy;                  T_rec_norm_o = 6.9+do;
T_icu_y = 2.5;                          T_icu_o = 1.5;          
T_rec_icu_y = 12.3-T_icu_y-dy;          T_rec_icu_o = 12.3-T_icu_o+do;    
T_vent_y = 2.5;                         T_vent_o = 1.5;         
T_rec_vent_y = 14.5-T_vent_y-T_icu_y-dy;T_rec_vent_o = 14.5-T_vent_o-T_icu_o+do;    
lambda_n_y = 2.32/100;                  lambda_n_o = 31.86/100;
lambda_c_y = 0.21;                      lambda_c_o = 0.31;
lambda_v_y = 0.5;                       lambda_v_o = 0.65;
lambda_d_y = 0.0232/
% infectious I
alpha_y = lambda_n_y./T_hosp_y;         alpha_o = lambda_n_o./T_hosp_o;          alpha = rho*alpha_o+(1-rho)*alpha_y;
alpha_r_y = (1-lambda_n_y)/T_inf_y;     alpha_r_o = (1-lambda_n_o)/T_inf_o;      zeta_i = rho*alpha_r_o+(1-rho)*alpha_r_y;
% normal hospital N
theta_n = rho/(1-rho)*lambda_n_o/lambda_n_y;    theta_n = theta_n/(1+theta_n);
beta_y = lambda_c_y/T_icu_y;            beta_o = lambda_c_o/T_icu_o;             beta = theta_n*beta_o+(1-theta_n)*beta_y;
beta_r_y = (1-lambda_c_y)/T_rec_norm_y; beta_r_o = lambda_c_o/T_rec_norm_o;      zeta_n = theta_n*beta_r_o+(1-theta_n)*beta_r_y;
% intensive care C
theta_c = theta_n*lambda_c_o/lambda_c_y;
gamma_y = lambda_v_y/T_vent_y;          gamma_o = lambda_v_o/T_vent_o;           gamma = theta_c*gamma_o+(1-theta_c)*gamma_y;
gamma_r_y = (1-lambda_v_y)/T_rec_icu_y; gamma_r_o = lambda_v_o/T_rec_icu_o;      zeta_c = theta_c*gamma_r_o+(1-theta_c)*gamma_r_y;
% ventilation V
theta_v = theta_c*lambda_v_o/lambda_v_y;
delta_y = lambda_d_y/T_death_y;         delta_o = lambda_d_o/T_death_o;          delta = theta_v*delta_o+(1-theta_v)*delta_y;
delta_r_y = (1-lambda_d_y)/T_rec_vent_y;delta_r_o = lambda_d_o/T_rec_vent_o;     zeta_v = theta_v*delta_r_o+(1-theta_v)*delta_r_y;

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

% initialization
dI_data = smooth_series(x.NewCases(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
D = smooth_series(x.Deaths(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
V = smooth_series(h.Ventilation(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
C = smooth_series(h.ICU(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
H = smooth_series(h.Hospitalizations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
N = H-C-V;

% calculation
V_D = smooth_series(D(2:end)-D(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
delta_t = V_D./V(1:end-1); k_v_t = delta_t./delta;
zeta_v_t = theta_v*(1-omega_y*k_v_t)/T_rec_vent_y+(1-theta_v)*(1-omega_o*k_v_t)/T_rec_vent_o;
V_R = zeta_v_t.*V(1:end-1);
C_V = V(2:end)-V(1:end-1)+V_D+V_R;
gamma_t = C_V./C(1:end-1); k_c_t = gamma_t./gamma;
zeta_c_t = theta_n*(1-omega_y*k_c_t)/T_rec_icu_y+(1-theta_n)*(1-omega_o*k_c_t)/T_rec_icu_o;
C_R = zeta_c_t.*C(1:end-1);
N_C = C(2:end)-C(1:end-1)+C_V+C_R;
beta_t = N_C./C(1:end-1); k_n_t = beta_t./beta;
zeta_n_t = theta_n*(1-omega_y*k_n_t)/T_rec_norm_y+(1-theta_n)*(1-omega_o*k_n_t)/T_rec_norm_o;
N_R = zeta_n_t.*N(1:end-1);
I_N = N(2:end)-N(1:end-1)+N_C+N_R;
I = I_N./alpha(1:T-1);
% I_R = alpha_r.*I;
X = I(2:end)-I(1:end-1).*(1-alpha_r-alpha(1:T-2));


% calculation
% H_D = smooth_series(D(2:end)-D(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
% beta_d_t = H_D./H(1:end-1); k_v_t = beta_d_t./beta_d;
% beta_r_t = theta*(1-omega_y*k_v_t)/T_rec_y+(1-theta)*(1-omega_o*k_v_t)/T_rec_o;
% % H_R = beta_r_t.*H(1:end-1);
% I_H = H(2:end)-H(1:end-1).*(1-beta_d_t-beta_r_t);
% I = I_H./alpha_h(1:T-1);
% % I_R = alpha_r.*I;
% X = I(2:end)-I(1:end-1).*(1-alpha_r-alpha_h(1:T-2));

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
p.T_delay = T_delay;
p.T_inf_y = T_inf_y;
p.T_inf_o = T_inf_o;
p.T_hosp_y = T_hosp_y;
p.T_hosp_o = T_hosp_o;
p.T_death_y = T_death_y; 
p.T_death_o = T_death_o; 
p.T_rec_y = T_rec_y; 
p.T_rec_o = T_rec_o;
p.alpha_h_y = alpha_y;
p.alpha_h_o = alpha_o;
p.alpha_r_y = alpha_r_y;
p.alpha_r_o = alpha_r_o;
p.lambda_n_y = lambda_n_y;
p.lambda_n_o = lambda_n_o;
p.beta_d_t = beta_d_t;
p.beta_r_t = beta_r_t;
p.k_d_t = k_v_t;
p.rho = rho;
p.theta = theta;
p.beta_d = beta_d;
p.beta_d_o = beta_d_o;
p.beta_d_y = beta_d_y;
p.omega_y = omega_y;
p.omega_o = omega_o;
% p.T_icu = T_icu; p.T_rec_icu = T_rec_icu;
% p.T_vent = T_vent; p.T_rec_vent = T_rec_vent;
% p.T_death = T_death; 
% p.alpha_d = alpha_d; p.alpha_v = alpha_v; p.alpha_c = alpha_c;
% p.lambda = lambda;

    function [x] = adjust_tail(x,k)
        dx = x(T-k)-x(T-k-1);
        x(T-k+1) = x(T-k)+2/3*dx;
        x(T-k+2) = x(T-k+1)+1/3*dx;
        for j=3:k
            x(T-k+j) = x(T-k+j-1)+1/3*1/(j-1)*dx;
        end        
    end

%     function [x] = get_rv(y)
%         shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
%         L = length(shape0);
%         shape0_vec = repmat(shape0,N,1);
%         scale0_vec = scale0*ones(N,L);
%         x = gamrnd(shape0_vec,scale0_vec);
%     end

end