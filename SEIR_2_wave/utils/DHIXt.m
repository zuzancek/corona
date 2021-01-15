function [X,I,obs_ratio_adj,sa,p] = DHIXt(x,h,d,s,dateFrom,dateTo,t0,t1,params,delay)

T = dateTo-dateFrom+1;
method_data = s.smoothing_method_data; 
method_params = s.smoothing_method_params;
firstData = dateFrom-31;% params.dataFrom;
tshift = dateFrom-firstData;
idxstr.D.validFrom = 1;
idxstr.D.dispFrom = tshift+1;
idxstr.H.validFrom = 1;
idxstr.H.dispFrom = tshift+1;
idxstr.HD.validFrom = 1;
idxstr.HD.dispFrom = tshift+1;

varsigma = extend(double(params.death_old_ratio),tshift);
% cfr_hospitals = method(params.cfr_hospitals);
% delta = cfr_hospitals(dateFrom:dateTo);
rho = method_params(params.cases_old_ratio(firstData:dateTo));
asymp_ratio = method_params(params.asymp_ratio(firstData:dateTo));
sigma = method_params(params.asymp_ratio(firstData:dateTo));

% delay in testing (gradual)
T_delay = NaN+zeros(T,1); T_delay(1) = 0;
dlen = length(delay.v);
for i=1:dlen
    T_delay(delay.at(i)-dateFrom) = delay.v(i);
end
if dlen
    T_delay(T) = delay.v(end);
end
T_delay = method_params(interp1(find(~isnan(T_delay)),T_delay(find(~isnan(T_delay))),1:T)'); %#ok<FNDSB>

% death: najpr time-inconsistent, potom konzistentne
omega_y = s.omega_y;       omega_o = s.omega_o;         omega = (omega_o.*varsigma+omega_y)./(1+varsigma);
T_death_y = s.T_death_y;   T_death_o = s.T_death_y;     T_death = (T_death_o.*varsigma+T_death_y)./(1+varsigma);
k_death = 10; x_death = 1:k_death; T_death_shape = T_death*s.T_death.std^2; T_death_scale = 1/s.T_death.std^2;
p_T_death = pdf('Gamma',repmat(x_death,length(varsigma),1),repmat(T_death_shape,1,k_death),repmat(T_death_scale,length(varsigma),k_death));
% hospitalizations
% old-young share in hospitals
zeta0 = (1-varsigma)./varsigma.*omega_y./omega_o; zeta = zeta0./(1+zeta0);
% recovery     
T_rec_y = s.T_rec_y;    T_rec_o = s.T_rec_o;            T_rec = zeta.*T_rec_o+(1-zeta).*T_rec_y;
k_rec = 20; x_rec = 1:k_rec;   T_rec_shape = T_rec*s.SI.std^2; T_rec_scale = 1/s.SI.std^2;
p_T_rec = pdf('Gamma',repmat(x_rec,length(zeta),1),repmat(T_rec_shape,1,k_rec),repmat(T_rec_scale,length(zeta),k_rec));

% T_death_y = s.T_death_y;                             alpha_hdy = omega_y/T_death_y; s.alpha_hdy = alpha_hdy;
% T_death_o = s.T_death_o;                             alpha_hdo = omega_o/T_death_o; s.alpha_hdo = alpha_hdo;
% T_rec_y = s.T_rec_y;                         alpha_hry = (1-omega_y)./T_rec_y;
% T_rec_o = s.T_rec_o;                         alpha_hro = (1-omega_o)./T_rec_o;
% eta_y = s.eta_y;                                     T_hosp_y = s.T_hosp_y;  alpha_ihy = eta_y./T_hosp_y;   
% eta_o = s.eta_o;                                     T_hosp_o = s.T_hosp_o;  alpha_iho = eta_o./T_hosp_o;  
% T_sick_y = s.T_sick_y.mean-s.T_test.mean-T_delay;    alpha_iry = (1-eta_y)./T_sick_y;
% T_sick_o = s.T_sick_o.mean-s.T_test.mean-T_delay;    alpha_iro = (1-eta_o)./T_sick_o;

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);     
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
% D(t+1) = D(t)+H_D(t);

% initialization
dI_data = method_data(x.NewCases(dateFrom:dateTo));
D = x.Deaths(firstData:dateTo)*d(dateFrom)/x.Deaths(dateFrom); 
D(tshift+1:end) = method_data(d(dateFrom:dateTo));
H = method_data(h.Hospitalizations(firstData:dateTo));

H0 = method_data(h.Hospitalizations(dateFrom-k_death:dateTo));
H01 = method_data(h.Hospitalizations(dateFrom-k_rec:dateTo));

HD = method_data(D(2:end)-D(1:end-1));  HD = [HD(1);HD(:)];
hd = get_wa(p_T_death(k_death+1:end,:),H,omega(k_death+1:end),k_death);
gamma_hd =  extend(HD(k_death+1:end)./hd,k_death);
gamma_hd = method_params(gamma_hd);
alpha_rec = 1-omega.*gamma_hd; 
HR = alpha_rec./T_rec.*H;
% HR = get_wa_inv(p_T_rec,H01,alpha_rec,k_rec);


% H_y_H_o = adjust_series(alpha_hdo./alpha_hdy.*H_D_y./H_D_o);
% H_y = method_params(H_y_H_o./(H_y_H_o+1)).*H;
% H_o = H-H_y;
% % alpha_hdy = method_params(H_D_y./H_y(1:end-1)); alpha_hry = method_params((1-alpha_hdy.*T_death_y)./T_rec_y);alpha_hry = [alpha_hry(1);alpha_hry];
% % alpha_hdo = method_params(H_D_o./H_y(1:end-1)); alpha_hro = method_params((1-alpha_hdo.*T_death_o)./T_rec_o);alpha_hro = [alpha_hro(1);alpha_hro];
% H_R_o = method_data(alpha_hro.*H_o);
% I_H_o = method_data(H_o(2:end)-H_o(1:end-1)+H_D_o+H_R_o(1:end-1));
% H_R_y = method_data(alpha_hry.*H_y);
% I_H_y = method_data(H_y(2:end)-H_y(1:end-1)+H_D_y+H_R_y(1:end-1));
% I_o = method_data(I_H_o./alpha_iho);
% % I_R_o = method_data(alpha_iro(1:end-1).*I_o); I_o = [I_o(1);I_o];
% X_o = method_data(I_o(2:end)-I_o(1:end-1))+I_H_o+I_R_o;
% I_y = method_data(I_H_y./alpha_ihy);
% I_R_y = method_data(alpha_iry(1:end-1).*I_y); I_y = [I_y(1);I_y];
% X_y = method_data(I_y(2:end)-I_y(1:end-1))+I_H_y+I_R_y;
% X = X_o+X_y;
% I = I_o+I_y;
rho_real = method_params(X_o./X);rho_real = [rho_real(1);rho_real];
plot(X,'linewidth',1);hold on;plot(dI_data,'k','linewidth',1);grid on;
figure;plot(rho);hold on;plot(rho_real);

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
idx = dateFrom:max(idx);
X(idx) = dI_data_reported(idx);
% X(idx) = dI_data_reported(idx); X(dateFrom:min(idx)) = dI_data_reported(dateFrom:min(idx));
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
p.rho = rho_real;
p.rho_ext = rho_ext;
p.varsigma = varsigma;

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
        ys = method_data(y);
    end

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

    function [x] = get_wa(weight,Z,alpha,idxFrom)
        sz = size(weight);
        k = sz(2);k0 = sz(1);
        %k = length(weight);
        t = length(Z)-idxFrom;
        if k0==1
            W = repmat(weight,t,1);
            A = repmat(alpha./(1:k),t,1);
        else
            W = weight;
            A = alpha./(1:k);
        end
        J = repmat(1:k,t,1)+repmat((0:t-1)',1,k);
        L = 0*(k-1)+repmat((1:t)',1,k);
        Weight_mat = sparse(L,J,W);
        Alpha_mat = sparse(L,J,A);
        Weight_mat = Weight_mat./sum(Weight_mat,2);
        x = (Weight_mat.*Alpha_mat)*Z(end-t-k+1:end-1);
    end

    function [x] = get_wa_inv(weight,Z,phi,idxFrom)
        sz = size(weight);
        k = sz(2);k0 = sz(1);
        t = length(Z)-idxFrom;
        a = phi./(1:k);
        if k0==1
            W = repmat(weight,t,1);
            A = repmat(a,t,1);
        else
            W = weight;
            A = a;
        end        
        J = repmat(1:k,t,1)+repmat((0:t-1)',1,k);
        L = (k-1)+repmat((1:t)',1,k);
        U0 = tril(repmat(1:k-1,k-1,1)); 
        U0(U0==0) = k+1; w(end+1) = 0;
        J0 = repmat(1:k-1,k-1,1);
        L0 = repmat((1:k-1)',1,k-1);
        W0 = w(U0(U0~=0));
        A0 = a(U0(U0~=0));
        Weight_mat = sparse([L(:),L0(:)],[J(:),J0(:)],[W(:),W0(:)]);
        Alpha_mat = sparse([L(:),L0(:)],[J(:),J0(:)],[A(:),A0(:)]);
        Weight_mat = Weight_mat./sum(Weight_mat,2);
        x = (Weight_mat.*Alpha_mat)\Z(end-t-k+1:end-1); 
    end

    function [y] = extend(x,t0)
        xlen = length(x);
        z = x(1)+zeros(xlen+t0,1);
        z(t0+1:end) = x;
        y = method_params(z);
    end

end