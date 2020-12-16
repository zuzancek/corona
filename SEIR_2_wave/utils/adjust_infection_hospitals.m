function [X,I,obs_ratio_adj,sa,p] = adjust_infection_hospitals(x,h,s,dateFrom,dateTo,t0,t1,sigma,omega,cfr,delay)

T = dateTo-dateFrom+1;

rho = s.old_share;

% delay in testing (gradual)
T_delay_0 = delay.v0;               T_delay_1 = delay.v1;               
T_delay_at = delay.at;
T_delay = zeros(T,1)+T_delay_0;     T_delay(T_delay_at-dateFrom:end) = T_delay_1;
T_delay = smooth_series(T_delay,s.smooth_width,s.smooth_type,s.smooth_ends);

T_inf = s.T_inf.mean;               T_test = (1+s.T_pre.mean)+T_delay;
T_inf_y = T_inf;                    T_inf_o = T_inf+1;
T_hosp_y = 7.02+T_test;             T_hosp_o = 3.24+T_test;         
lambda_y = 2.32/100;                lambda_o = 28.86/100;
alpha_h_y = lambda_y./T_hosp_y;     alpha_h_o = lambda_o./T_hosp_o;        alpha_h = rho*alpha_h_y+(1-rho)*alpha_h_o;
alpha_r_y = (1-lambda_y)/T_inf_y;   alpha_r_o = (1-lambda_o)/T_inf_o;      alpha_r = rho*alpha_r_o+(1-rho)*alpha_r_y;
T_death_y = 3.41;                   T_death_o = 4.59;
T_rec_y = 4.56;                     T_rec_o = 5.65;
omega_y = 5.15/100;                 omega_o = 37.16;
theta = rho/(1-rho)*lambda_o/lambda_y; 
theta = theta/(1+theta);
beta_d_y = omega_y/T_death_y;       beta_d_o = omega_o/T_death_o;          beta_d = theta*beta_d_o+(1-theta)*beta_d_y;

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);     
%       I_H(t) = alpha_h*I(t);   I_R(t) = alpha_r*I(t);
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
%       H_D(t) = beta_d*H(t);    H_R(t) = beta_r*N(t);
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
p.alpha_h_y = alpha_h_y;
p.alpha_h_o = alpha_h_o;
p.alpha_r_y = alpha_r_y;
p.alpha_r_o = alpha_r_o;
p.lambda_y = lambda_y;
p.lambda_o = lambda_o;
p.beta_d_t = beta_d_t;
p.beta_r_t = beta_r_t;
p.k_d_t = k_d_t;
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