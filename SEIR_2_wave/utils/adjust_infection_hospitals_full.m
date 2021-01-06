function [X,I,obs_ratio_adj,sa,p] = adjust_infection_hospitals_full(x,h,d,s,dateFrom,dateTo,t0,t1,params,delay,srec)

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
T_short_0 = srec.v0;                T_short_1 = 0*srec.v1;           
T_short_at = srec.at;
T_short = zeros(T,1)+T_short_0;     T_short(T_short_at-dateFrom:end) = T_short_1;
T_short = method_params(T_short);

% parameters
% I->H/N
lambda_iny = s.eta_y+0.04;    T_hosp_y = s.T_hosp_y;   T_sick_y = s.T_sick_y.mean-s.T_test.mean-T_delay;    
lambda_ino = s.eta_o+0.04;    T_hosp_o = s.T_hosp_o;   T_sick_o = s.T_sick_o.mean-s.T_test.mean-T_delay;    
alpha_iry = (1-lambda_iny)./T_sick_y; alpha_iny = lambda_iny./T_hosp_y; 
alpha_iro = (1-lambda_ino)./T_sick_o; alpha_ino = lambda_ino./T_hosp_o; 
% N->C
lambda_ncy = (7.61-3)*3/100;    T_icu_y = 3;      T_rec_ny = 6.9-2-T_short;
lambda_nco = (7.61+3)*3/100;    T_icu_o = 3;      T_rec_no = 6.9+2-T_short;
alpha_ncy = lambda_ncy./T_icu_y; alpha_nry = (1-lambda_ncy)./T_rec_ny; 
alpha_nco = lambda_nco./T_icu_o; alpha_nro = (1-lambda_nco)./T_rec_no; 
% C->V
lambda_cvy = (9.86-3)*3/100;    T_vent_y = 3;      T_rec_cy = 12.3-2-T_icu_y;
lambda_cvo = (9.86+3)*3/100;    T_vent_o = 3;      T_rec_co = 12.3+3-T_icu_o;
alpha_cvy = lambda_cvy./T_vent_y; alpha_cry = (1-lambda_cvy)./T_rec_cy; 
alpha_cvo = lambda_cvo./T_vent_o; alpha_cro = (1-lambda_cvo)./T_rec_co; 
% V->D
T_death_y = 8.5-T_vent_y-T_icu_y;     lambda_vdy = (9.73-3)*T_death_y/100;     T_rec_cy = 15.8-2-T_icu_y-T_vent_y;
T_death_o = 12.6-T_vent_o-T_icu_o;    lambda_vdo = (9.73+3)*T_death_o/100;     T_rec_co = 15.8+3-T_icu_o-T_vent_o;
alpha_vdy = lambda_vdy./T_death_y;    alpha_vry = (1-lambda_cvy)./T_rec_cy; 
alpha_vdo = lambda_vdo./T_death_o;    alpha_vro = (1-lambda_cvo)./T_rec_co; 

%%
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

%%
% initialization
dI_data = method_data(x.NewCases(dateFrom:dateTo));
D = method_data(d(dateFrom:dateTo));
H = method_data(h.Hospitalizations(dateFrom:dateTo));
V = method_data(h.Ventilation(dateFrom:dateTo));
C = method_data(h.ICU(dateFrom:dateTo));
N = H-C-V;

% calculation
V_D = method_data(D(2:end)-D(1:end-1));
V_D_o = varsigma(1:end-1).*V_D;
V_D_y = V_D-V_D_o;
%
V_y_V_o = adjust_series(alpha_vdo./alpha_vdy.*V_D_y./V_D_o);
V_y = V_y_V_o./(V_y_V_o+1).*V;
V_o = V-V_y;
V_R_o = alpha_vro.*V_o;
V_R_y = alpha_vry.*V_y;
C_V_o = method_data(V_o(2:end)-V_o(1:end-1))+V_D_o+V_R_o(1:end-1);
C_V_y = method_data(V_y(2:end)-V_y(1:end-1))+V_D_y+V_R_y(1:end-1);
%
C_y_C_o = adjust_series(alpha_cvo./alpha_cvy.*C_V_y./C_V_o);
C_y = C_y_C_o./(C_y_C_o+1).*C;
C_o = C-C_y;
C_R_o = alpha_cro.*C_o;
C_R_y = alpha_cry.*C_y;
N_C_o = method_data(C_o(2:end)-C_o(1:end-1))+C_V_o+C_R_o(1:end-1);
N_C_y = method_data(C_y(2:end)-C_y(1:end-1))+C_V_y+C_R_y(1:end-1);
%
N_y_N_o = adjust_series(alpha_nco./alpha_ncy.*N_C_y./N_C_o);
N_y = N_y_N_o./(N_y_N_o+1).*N;
N_o = N-N_y;
N_R_o = alpha_nro.*N_o;
N_R_y = alpha_nry.*N_y;
I_N_o = method_data(N_o(2:end)-N_o(1:end-1))+N_C_o+N_R_o(1:end-1);
I_N_y = method_data(N_y(2:end)-N_y(1:end-1))+N_C_y+N_R_y(1:end-1);
%
I_o = I_N_o./alpha_ino;
I_y = I_N_y./alpha_iny;
I_R_o = alpha_iro(1:end-1).*I_o; I_o = [I_o(1);I_o];
I_R_y = alpha_iry(1:end-1).*I_y; I_y = [I_y(1);I_y];
X_o = method_data(I_o(2:end)-I_o(1:end-1))+I_N_o+I_R_o;
X_y = method_data(I_y(2:end)-I_y(1:end-1))+I_N_y+I_R_y;
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
p = struct();
p.alpha_iry = alpha_iry; p.alpha_iny = alpha_iny; 
p.alpha_iro = alpha_iro; p.alpha_ino = alpha_ino; 
p.alpha_ncy = alpha_ncy; p.alpha_nry = alpha_nry; 
p.alpha_nco = alpha_nco; p.alpha_nro = alpha_nro; 
p.alpha_cvy = alpha_cvy; p.alpha_cry = alpha_cry; 
p.alpha_cvo = alpha_cvo; p.alpha_cro = alpha_cro; 
p.alpha_vdy = alpha_vdy; p.alpha_vry = alpha_vry; 
p.alpha_vdo = alpha_vdo; p.alpha_vro = alpha_vro; 
p.lambda_iny = lambda_iny; p.lambda_ino = lambda_ino;
p.lambda_ncy = lambda_ncy; p.lambda_nco = lambda_nco;
p.lambda_cvy = lambda_cvy; p.lambda_cvo = lambda_cvo;

    function [x] = adjust_tail(x,k)
        dx = x(T-k)-x(T-k-1);
        x(T-k+1) = x(T-k)+2/3*dx;
        x(T-k+2) = x(T-k+1)+1/3*dx;
        for j=3:k
            x(T-k+j) = x(T-k+j-1)+1/3*1/(j-1)*dx;
        end        
    end

    function [y,ys] = adjust_series(x)
        y = NaN+zeros(T,1);
        y(1:length(x)) = x;
        if isnan(y(1))
            y(1) = mean(x,'omitnan');
        end
        if isnan(y(end))
            if isnan(x(end))
                y(end) = mean(x,'omitnan');
            else
                y(end) = x(end);
            end
        end
        idx_y = ~isnan(y);
        y = interp1(find(idx_y),y(idx_y),1:T,'spline')';
        ys = smooth_series(y);
    end

%     function [x] = get_rv(y)
%         shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
%         L = length(shape0);
%         shape0_vec = repmat(shape0,N,1);
%         scale0_vec = scale0*ones(N,L);
%         x = gamrnd(shape0_vec,scale0_vec);
%     end

end