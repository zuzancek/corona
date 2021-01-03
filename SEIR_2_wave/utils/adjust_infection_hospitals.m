function [X,I,obs_ratio_adj,sa,p] = adjust_infection_hospitals(x,h,d,s,dateFrom,dateTo,t0,t1,params,delay,srec)

T = dateTo-dateFrom+1;
method_data = s.smoothing_method_data; 
method_params = s.smoothing_method_params; 

death_old_ratio = method_params(params.death_old_ratio);
varsigma = death_old_ratio(dateFrom:dateTo);
% cfr_hospitals = method(params.cfr_hospitals);
% delta = cfr_hospitals(dateFrom:dateTo);
cases_old_ratio = method_params(params.cases_old_ratio);
rho = cases_old_ratio(dateFrom:dateTo);
rho_ext = cases_old_ratio;
asymp_ratio = method_params(params.asymp_ratio);
sigma = asymp_ratio(dateFrom:dateTo);

% delay in testing (gradual)
T_delay_0 = delay.v0;               T_delay_1 = delay.v1;               
T_delay_at = delay.at;
T_delay = zeros(T,1)+T_delay_0;     T_delay(T_delay_at-dateFrom:end) = T_delay_1;
T_delay = method_params(T_delay);
% shortening recovery peiod
T_short_0 = srec.v0;                T_short_1 = srec.v1;           
T_short_at = srec.at;
T_short = zeros(T,1)+T_short_0;     T_short(T_short_at-dateFrom:end) = T_short_1;
T_short = method_params(T_short);

omega_y = s.omega_y;
omega_o = s.omega_o;
T_death_y = s.T_death_y;        alpha_hdy = omega_y/T_death_y;
T_death_o = s.T_death_o;        alpha_hdo = omega_o/T_death_o;
T_rec_y = s.T_rec_y-T_short;    alpha_hry = (1-omega_y)./T_rec_y;
T_rec_o = s.T_rec_o-T_short;    alpha_hro = (1-omega_o)./T_rec_o;
eta_y = s.eta_y;                alpha_ihy = eta_y;            
T_hosp_y = s.T_hosp_y;
eta_o = s.eta_o;                alpha_iho = eta_o;            
T_hosp_o = s.T_hosp_o;
T_sick_y = s.T_sick_y.mean-s.T_test.mean-T_delay;         alpha_iry = (1-eta_y*T_hosp_y)./T_sick_y;
T_sick_o = s.T_sick_o.mean-s.T_test.mean-T_delay;         alpha_iro = (1-eta_o*T_hosp_o)./T_sick_o;

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);     
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
% D(t+1) = D(t)+H_D(t);

% initialization
dI_data = method_data(x.NewCases(dateFrom:dateTo));
D = method_data(d(dateFrom:dateTo));
H = method_data(h.Hospitalizations(dateFrom:dateTo));

H_D = method_data(D(2:end)-D(1:end-1));
H_D_o = varsigma(1:end-1).*H_D;
H_D_y = H_D-H_D_o;
H_o = H_D_o./alpha_hdo; H_o(end+1) = H_o(end);
H_y = H-H_o;
H_R_o = alpha_hro(1:end-1).*H_o(1:end-1);
I_H_o = method_data(H_o(2:end)-H_o(1:end-1))+H_D_o+H_R_o;
I_o = I_H_o./alpha_iho;
I_R_o = alpha_iro(1:end-1).*I_o; I_o = [I_o(1);I_o];
X_o = method_data(I_o(2:end)-I_o(1:end-1))+I_H_o+I_R_o;
H_R_y = alpha_hry(1:end-1).*H_y(1:end-1);
I_H_y = method_data(H_y(2:end)-H_y(1:end-1))+H_D_y+H_R_y;
I_y = I_H_y./alpha_ihy;
I_R_y = alpha_iry(1:end-1).*I_y; I_y = [I_y(1);I_y];
X_y = method_data(I_y(2:end)-I_y(1:end-1))+I_H_y+I_R_y;
X = X_o+X_y;
I = I_o+I_y;

% adjust series endpoints and get ratio
obs_ratio_adj = tseries(t0:t1,s.obs_ratio);
X = adjust_tail(X,2);
X_o = adjust_tail(X_o,2);
X_y = adjust_tail(X_y,2);
X = tseries(dateFrom:dateTo,method_data(X));
X_o = tseries(dateFrom:dateTo,method_data(X_o));
X_y = tseries(dateFrom:dateTo,method_data(X_y));
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
sa.Xs = (1-sigma(1)).*X;
sa.Xo = X_o;
sa.Xy = X_y;
sa.Xa = X-sa.Xs;
sa.dIa_data_reported = dI_data_reported.*sigma;
sa.dIs_data_reported = dI_data_reported-sa.dIa_data_reported;
sa.loss_a = sa.Xa-sa.dIa_data_reported;
sa.loss_s = sa.Xs-sa.dIs_data_reported; idx = find(sa.loss_s<0); sa.loss_s(idx) = 0; %#ok<FNDSB>
sa.loss_o = sa.Xo-dI_data_reported_old;
sa.loss_y = sa.Xy-dI_data_reported_young;

% store params
p.alpha_hdy = alpha_hdy;
p.alpha_hdo = alpha_hdo;
p.alpha_hry = alpha_hry;
p.alpha_hro = alpha_hro;
p.alpha_ihy = alpha_ihy;
p.alpha_iho = alpha_iho;
p.alpha_iry = alpha_iry;
p.alpha_iro = alpha_iro;
p.T_rec_y = T_rec_y;
p.T_rec_o = T_rec_o;
p.T_sick_y = T_sick_y;
p.T_sick_o = T_sick_o;
p.rho = rho;
p.rho_ext = rho_ext;
p.varsigma = varsigma;

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