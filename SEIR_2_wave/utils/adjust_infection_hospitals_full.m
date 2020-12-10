function [X,I,obs_ratio_adj,sa] = adjust_infection_hospitals_full(x,h,s,dateFrom,dateTo,t0,t1,sigma)

T = dateTo-dateFrom+1;
T_hosp_data.init = 0; T_hosp_data.final = 0.5; T_hosp_data.bp = dd(2020,10,15);

%
T_symp = s.T_hosp.mean; 
T_hosp0 = T_hosp_data.init+T_symp; 
T_hosp1 = T_hosp_data.final+T_symp; 
T_hosp = T_hosp0+zeros(T-1,1);
bp = T_hosp_data.bp-dateFrom;
T_hosp(end) = T_hosp1; T_hosp(bp:end) = linspace(T_hosp0,T_hosp1,T-bp);
T_hosp = smooth_series(T_hosp,s.smooth_width,s.smooth_type,s.smooth_ends);

T_rec_norm = 6.9;
T_icu = 8-T_hosp;
T_rec_icu = 12.3-T_icu;
T_vent = 5-T_icu;
T_rec_vent = 14.5-(T_vent+T_icu);
T_death = 12-(T_vent+T_icu);
s.lambda = 6.78/100;
lambda = s.lambda; % 6.37/100;
alpha_in = lambda./T_hosp.*ones(T-1,1);
T_inf = s.SI.mean;
alpha_ir = (1-lambda)/T_inf;

% alpha_vd = 7.62/100;
% alpha_cv = 9.86/100;
% alpha_nc = 7.61/100;

dI_data = smooth_series(x.NewCases(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);

% definitions
D = smooth_series(x.Deaths(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
V = smooth_series(h.Ventilation(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
C = smooth_series(h.ICU(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
H = smooth_series(h.Hospitalizations(dateFrom:dateTo),s.smooth_width,s.smooth_type,s.smooth_ends);
N = H-C-V;

d_V_D = smooth_series(D(2:end)-D(1:end-1),s.smooth_width,s.smooth_type,s.smooth_ends);
alpha_vd = d_V_D./V(1:end-1);
alpha_d = T_death.*alpha_vd;
alpha_vr = (1-alpha_d)./T_rec_vent;
d_V_R = alpha_vr.*V(1:end-1);
d_V = V(2:end)-V(1:end-1);
d_C_V = d_V+d_V_D+d_V_R;
alpha_cv = d_C_V./C(1:end-1);
alpha_v = alpha_cv.*T_vent;
alpha_cr = (1-alpha_v)./T_rec_icu;
d_C_R = alpha_cr.*C(1:end-1);
d_C = C(2:end)-C(1:end-1);
d_N_C = d_C+d_C_V+d_C_R;
alpha_nc = d_N_C./N(1:end-1);
alpha_c = T_icu.*alpha_nc;
alpha_nr = (1-alpha_c)./T_rec_norm;
d_N_R = alpha_nr.*N(1:end-1);
d_N = N(2:end)-N(1:end-1);
d_I_N = d_N+d_N_C+d_N_R;
I = d_I_N./alpha_in;
d_I_R = I.*alpha_ir;
d_I = I(2:end)-I(1:end-1);
X = d_I+d_I_N(1:end-1)+d_I_R(1:end-1);

% adjust series endpoints and get ratio
obs_ratio_adj = tseries(t0:t1,s.obs_ratio);
X = adjust_tail(X,3);
X = tseries(dateFrom:dateTo,smooth_series(X,s.smooth_width,s.smooth_type,s.smooth_ends));
dI_data_real = resize(X,dateFrom:dateTo);
dI_data_reported = tseries(dateFrom:dateTo,dI_data);
delta = dI_data_reported./dI_data_real;

idx = find(dI_data_real<s.cases_min & dI_data_reported<s.cases_min & delta<1-s.ratio_threshold);
X(idx) = dI_data_reported(idx); X(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
dI_data_real(idx) = dI_data_reported(idx);dI_data_real(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
delta = dI_data_reported./dI_data_real;

obs_ratio_adj(dateFrom:dateTo) = smooth_series(delta*s.obs_ratio,s.smooth_width,s.smooth_type,s.smooth_ends);
XX = x.NewCases;
XX(dateFrom:dateTo) = X;
X = smooth_series(XX,s.smooth_width,s.smooth_type,s.smooth_ends);
sa = struct;
sa.Xs = (1-sigma(s.wave_2_from)).*X;
sa.Xa = X-sa.Xs;
sa.dIa_data_reported = dI_data_reported.*sigma(dateFrom:dateTo);
sa.dIs_data_reported = dI_data_reported-sa.dIa_data_reported;
sa.loss_a = sa.Xa-sa.dIa_data_reported;
sa.loss_s = sa.Xs-sa.dIs_data_reported;

    function [x] = adjust_tail(x,k)
        dx = x(T-k)-x(T-k-1);
        x(T-k+1) = x(T-k)+2/3*dx;
        x(T-k+2) = x(T-k+1)+1/3*dx;
        for j=3:k
            x(T-k+j) = x(T-k+j-1)+1/3*1/(j-1)*dx;
        end        
    end

end