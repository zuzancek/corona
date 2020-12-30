function [X,I,obs_ratio_adj,sa,p] = adjust_infection_hospitals_det(x,h,d,s,dateFrom,dateTo,t0,t1,sigma,omega,mort,delay)
% (x,h,d,s,dateFrom,dateTo,t0,t1,sigma,omega,cfr,delay)

T = dateTo-dateFrom+1;
method = @smooth_series; 
d_share_old = mort.share_old;
cfr = mort.cfr_hosp;
% method = s.smoothing_method;

rho = method(omega(dateFrom:dateTo)); %s.old_share;
rho_ext = method(omega); %s.old_share;
p = struct();

% delay in testing (gradual)
T_delay_0 = delay.v0;               T_delay_1 = delay.v1;               
T_delay_at = delay.at;
T_delay = zeros(T,1)+T_delay_0;     T_delay(T_delay_at-dateFrom:end) = T_delay_1;
T_delay = method(T_delay);

omega_y = 5.15/100;         T_death_y = 3.41;      alpha_hdy = omega_y/(T_death_y);
omega_o = (37.16/100);      T_death_o = 4.59;      alpha_hdo = omega_o/(T_death_o);
                            T_rec_y = 5;           alpha_hry = (1-omega_y)./T_rec_y;
                            T_rec_o = 7;           alpha_hro = (1-omega_o)./T_rec_o;
kappa = 1;
eta_y = 2.32/100;           alpha_ihy = kappa*eta_y;
eta_o = 31.86/100;          alpha_iho = kappa*eta_o;
T_inf_y = s.SI.mean;        alpha_iry = (1-kappa*eta_y)./T_inf_y;
T_inf_o = T_inf_y+1;        alpha_iro = (1-kappa*eta_o)./T_inf_o;

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);     
%       I_H(t) = alpha_h*I(t);   I_R(t) = alpha_r*I(t);
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
%       H_D(t) = beta_d*H(t);    H_R(t) = beta_r*N(t);
% D(t+1) = D(t)+H_D(t);

% initialization
% s.smooth_width = 7;
dI_data = method(x.NewCases(dateFrom:dateTo));
D = method(d(dateFrom:dateTo));%smooth_series(x.Deaths(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
H = method(h.Hospitalizations(dateFrom:dateTo));

H_D = method(D(2:end)-D(1:end-1));
alpha_hd = method(H_D./H(1:end-1));
theta = (alpha_hd-alpha_hdy)./(alpha_hdo-alpha_hdy);
H_D_o = theta.*H_D;
H_o = H_D_o./alpha_hdo; H_o(end+1) = H_o(end);
H_R_o = alpha_hro.*H_o(1:end-1);
I_H_o = method(H_o(2:end)-H_o(1:end-1))+H_D_o+H_R_o;
I_o = I_H_o./alpha_iho;
I_R_o = alpha_iro.*I_o; I_o = [I_o(1);I_o];
X_o = method(I_o(2:end)-I_o(1:end-1))+I_H_o+I_R_o;
H_D_y = H_D-H_D_o;
H_y = H-H_o;
H_R_y = alpha_hry.*H_y(1:end-1);
I_H_y = method(H_y(2:end)-H_y(1:end-1))+H_D_y+H_R_y;
I_y = I_H_y./alpha_ihy;
I_R_y = alpha_iry.*I_y; I_y = [I_y(1);I_y];
X_y = method(I_y(2:end)-I_y(1:end-1))+I_H_y+I_R_y;
X = X_o+X_y;

% adjust series endpoints and get ratio
obs_ratio_adj = tseries(t0:t1,s.obs_ratio);
X = adjust_tail(X,2);
X = tseries(dateFrom:dateTo,method(X));
%X_o = adjust_tail(X_o,2);
X_o = tseries(dateFrom:dateTo,method(X_o));
X_y = adjust_tail(X_y,2);
X_y = tseries(dateFrom:dateTo,method(X_y));
dI_data_real = resize(X,dateFrom:dateTo);
dI_data_reported = tseries(dateFrom:dateTo,dI_data);
dI_data_reported_old = dI_data_reported.*rho;
dI_data_reported_young = dI_data_reported-dI_data_reported_old;
delta = dI_data_reported./dI_data_real;

idx = find(dI_data_real<s.cases_min & dI_data_reported<s.cases_min & delta<1-s.ratio_threshold);
X(idx) = dI_data_reported(idx); X(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
X_o(idx) = dI_data_reported_old(idx); X_o(dateFrom:min(idx)) = dI_data_reported_old(dateFrom:min(idx));
X_y(idx) = dI_data_reported_old(idx); X_y(dateFrom:min(idx)) = dI_data_reported_young(dateFrom:min(idx));
dI_data_real(idx) = dI_data_reported(idx);dI_data_real(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
delta = dI_data_reported./dI_data_real;

obs_ratio_adj(dateFrom:dateTo) = smooth_series(delta*s.obs_ratio,s.smooth_width,s.smooth_type,s.smooth_ends);
XX = x.NewCases;
XX(dateFrom:dateTo) = X;
X = smooth_series(XX,s.smooth_width,s.smooth_type,s.smooth_ends);
XX = x.NewCases.*rho_ext;
XX(dateFrom:dateTo) = X_o;
X_o = smooth_series(XX,s.smooth_width,s.smooth_type,s.smooth_ends);
XX = x.NewCases.*(1-rho_ext);
XX(dateFrom:dateTo) = X_y;
X_y = smooth_series(XX,s.smooth_width,s.smooth_type,s.smooth_ends);

sa = struct;
sa.Xs = (1-sigma(dateFrom)).*X;%
sa.Xo = X_o;
sa.Xy = X_y;
sa.Xa = X-sa.Xs;
sa.dIa_data_reported = dI_data_reported.*sigma(dateFrom:dateTo);
sa.dIs_data_reported = dI_data_reported-sa.dIa_data_reported;
sa.loss_a = sa.Xa-sa.dIa_data_reported;
sa.loss_s = sa.Xs-sa.dIs_data_reported; idx = find(sa.loss_s<0); sa.loss_s(idx) = 0; %#ok<FNDSB>
sa.loss_o = sa.Xo-dI_data_reported_old;
sa.loss_y = sa.Xy-dI_data_reported_young;


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