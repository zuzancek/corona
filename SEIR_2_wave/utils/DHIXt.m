function [X,I,obs_ratio_adj,sa,p] = DHIXt(x,h,d,s,dateFrom,dateTo,t0,~,params,delay)

%% initialization
T = dateTo-dateFrom+1;
method_data = s.smoothing_method_data; 
method_params = s.smoothing_method_params;
tshift = s.firstData_offset;
firstData = -tshift+dateFrom;
cut = params.cutoff;
T_test_to_result = 1;
adj = params.adj;

varsigma = extend(double(resize(params.death_old_ratio,dateFrom:dateTo)),tshift);
rho = method_params(params.cases_old_ratio(firstData:dateTo));
sigma = method_params(params.asymp_ratio(firstData:dateTo));

%% testing
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

% death: najprv time-inconsistent
omega_y = s.omega_y;       omega_o = s.omega_o;          omega = (omega_o.*varsigma+omega_y)./(1+varsigma); 
T_death_y = s.T_death_y; % 1+1/0.1629;        
T_death_o = s.T_death_o; % 1+1/0.1092;
k_death = s.k_death; x_death = 1:k_death;
po = 1./(T_death_o.*varsigma); py = 1./(T_death_y.*(1-varsigma));
p_T_death = repmat(po.*py./(po-py),1,k_death).*(exp(-py*x_death)-exp(-po*x_death));
p_T_death = p_T_death./sum(p_T_death,2);
% hospitalizations
% old-young share in hospitals
zeta0 = (varsigma)./(1-varsigma).*omega_y./omega_o; zeta = zeta0./(1+zeta0);
% recovery in hospitals
T_rec_y = 7.9625;       T_rec_o = 11.7325;              T_rec = zeta.*T_rec_o+(1-zeta).*T_rec_y;
T_rec_std = 0.5547;
k_rec = s.k_rec;        x_rec = 1:k_rec;                T_rec_shape = T_rec*T_rec_std^2; T_rec_scale = 1/T_rec_std^2;
p_T_rec = pdf(s.T_rec_pdf_type,repmat(x_rec,length(zeta),1),repmat(T_rec_shape,1,k_rec),repmat(T_rec_scale,length(zeta),k_rec));
p_T_rec = p_T_rec./sum(p_T_rec,2);
% hospital admission
lambda_y = s.eta_y;    lambda_o = s.eta_o;      
theta0 = zeta0.*lambda_y./lambda_o; theta = theta0./(1+theta0); lambda = theta.*lambda_o+(1-theta).*lambda_y; 
k_hosp = s.k_hosp;     x_hosp = 1:k_hosp;  
py = 0.9425;           
pdf_T_hosp_y = -1./log(1-py).*py.^x_hosp./x_hosp; 
pdf_T_hosp_y = pdf_T_hosp_y./sum(pdf_T_hosp_y,2);
po = 0.8825;           pdf_T_hosp_o = -1./log(1-po).*po.^x_hosp./x_hosp; pdf_T_hosp_o = pdf_T_hosp_o./sum(pdf_T_hosp_o,2);
pdf_T_hosp = theta.*pdf_T_hosp_o+(1-theta).*pdf_T_hosp_y;  p_T_hosp = pdf_T_hosp./sum(pdf_T_hosp,2);
% T_hosp_y = s.T_hosp_y;        
% T_hosp_o = s.T_hosp_o;                   
% po = 1./(T_hosp_o(1).*theta); py = 1./(T_hosp_y(1).*(1-theta));
% p_T_hosp = repmat(po.*py./(po-py),1,k_hosp).*(exp(-py*x_hosp)-exp(-po*x_hosp));
% p_T_hosp = p_T_hosp./sum(p_T_hosp,2);
% infections
% recovery from sickness
T_delay = extend(T_delay,length(theta)-length(T_delay));
T_sick_y = s.T_sick_y;  T_sick_o = s.T_sick_o;   T_sick_std = s.T_sick_std; 
T_sick = theta.*T_sick_o+(1-theta).*T_sick_y-T_delay-0*T_test_to_result;
k_sick = s.k_sick;      x_sick = 1:k_sick;       T_sick_shape = T_sick*T_sick_std^2; T_sick_scale = 1/T_sick_std^2;
eta = 1-lambda;
p_T_sick = pdf(s.T_sick_pdf_type,repmat(x_sick,length(zeta),1),repmat(T_sick_shape,1,k_sick),repmat(T_sick_scale,length(zeta),k_sick));
p_T_sick = p_T_sick./sum(p_T_sick,2);

%% time shift (death->illness)
ks = s.t_shift_clin;        xs = 1:ks;
pp(1) = 1./(mean(T_death_o).*mean(varsigma));pp(2) = 1./(mean(T_death_y).*mean(1-varsigma));
% pp(3) = 1./(mean(T_hosp_o).*mean(theta));pp(4) = 1./(mean(T_hosp_y).*mean(1-theta));
lmat = repmat(pp,length(pp),1)-pp'; lmat(lmat==0) = NaN;
p_T_shift = prod(pp).*sum(exp(-pp'.*xs)./repmat(prod(lmat,2,'omitnan'),1,ks),1);
T_shift_death = ceil(dot(p_T_shift,xs));      % mean shift
[~,idx] = (max(theta));
T_shift_hosp = ceil(dot(xs(1:k_hosp)',p_T_hosp(idx,:)));
T_shift_sick = ceil(dot(xs(1:k_sick)',p_T_sick(idx,:)));
T_shift = T_shift_hosp+0*T_shift_death;
% ******* Equations
% I(t) = I(t-1)+X(t)-I_H(t)-I_R(t);     
% H(t) = H(t-1)+I_H(t)-H_D(t)-H_R(t);
% D(t) = D(t-1)+H_D(t);

% initialization
dI_data = method_data(x.NewCases(dateFrom:dateTo-cut));
dI_data_all = method_data(x.NewCases(dateFrom:dateTo));
D = x.Deaths(firstData:dateTo)*d(dateFrom)/x.Deaths(dateFrom); 
D(tshift+1:end) = method_data(d(dateFrom:dateTo));
H = method_data(h.Hospitalizations(firstData:dateTo));
AC = method_data(x.ActiveCases(firstData-k_hosp+2:dateTo));

% calculation
HD = method_data(D(2:end)-D(1:end-1));  
hd = method_params(get_wa(p_T_death,H,omega,k_death));
gamma_hd =  method_params(extend(HD(k_death+1:end)./hd(1:end-1),k_death));
omega = extend(method_params(omega(1:end-1).*gamma_hd),1);
HR = method_data(extend(get_wa(p_T_rec,H,1-omega,k_rec),k_rec));
IH = H(2:end)-H(1:end-1)+HR(1:end-1)+HD;
I = method_data(get_wa_inv(p_T_hosp(1:end-1,:),IH,AC,lambda(1:end-1),k_hosp));
IR = method_data(extend(get_wa(p_T_sick(1+0*T_shift_hosp:end-2,:),I,eta(1+0*T_shift_hosp:end-2),k_sick),k_sick));
dt = length(I(1:end-1-T_shift_sick));
X = method_data(I(2:end-T_shift_sick)-I(1:end-1-T_shift_sick)+IR(end-dt+1:end)...
    +IH(end-dt+1-(T_shift_sick-T_shift_hosp):end-(T_shift_sick-T_shift_hosp)));
X = extend_tail(X,3);

% X_orig = X;
kk = (1+adj/T_shift).^(0:T_shift-1); 
X(end-T_shift+1:end) = X(end-T_shift+1:end).*kk';

Xts = smooth_series(X(tshift:end)); Xts = tseries(dateFrom:dateFrom+length(Xts)-1,Xts);
Xrts = (X(tshift:end)); Xrts = tseries(dateFrom:dateFrom+length(Xrts)-1,Xrts);
Orts = tseries(dateFrom:dateFrom+length(dI_data)-1,dI_data);
Ots = smooth_series(Orts);
Orts0 = tseries(dateFrom:dateFrom+length(dI_data_all)-1,dI_data_all);
Ots0 = smooth_series(Orts0);

figure;
bar(Orts);hold on;
bar(Xrts);
bar(mov_median(params.h));

Len = length(Xts);
p = struct();
p.X_smooth = resize(Xts,dateFrom:dateFrom+Len-T_shift);
p.X_forecast_smooth = resize(Xts,dateFrom+Len-T_shift:dateFrom+Len-1);
p.X_raw = resize(Xrts,dateFrom:dateFrom+Len-T_shift);
p.X_forecast_raw = resize(Xrts,dateFrom+Len-T_shift:dateFrom+Len-1);
p.X_rep_smooth = Ots;
p.X_rep_forecast_smooth = resize(Ots0,enddate(Ots)+1:dateTo);
p.X_rep_raw = Orts;
rho_real_0 = method_params(lambda_y./lambda_o.*omega_y./omega_o.*varsigma./(1-varsigma));
rho_real_0 = [rho_real_0(1);rho_real_0];
rho_real = rho_real_0./(1+rho_real_0);

% adjust series endpoints and get ratio
X = X(tshift:end);
rho_real = rho_real(tshift+1:end);
sigma = sigma(tshift+1:end-cut);
dateTo_X = dateFrom+length(X)-1;
dateTo_R = dateFrom+length(dI_data)-1;
dateTo_0 = min(dateTo_X,dateTo_R);
obs_ratio_adj = tseries(t0:dateTo_0,s.obs_ratio);
X = Xts; %tseries(dateFrom:dateTo_X,method_data(X));
X_o = method_data(X.*rho_real);
X_y = X-X_o;
dI_data_real = resize(X,dateFrom:dateTo_X);
dI_data_reported = Ots;
dI_data_reported_old = dI_data_reported.*rho(tshift:tshift+length(dI_data_reported)-1);
dI_data_reported_young = dI_data_reported-dI_data_reported_old;
% delta = method_params(resize(dI_data_reported,dateFrom:dateTo_0)./resize(dI_data_real,dateFrom:dateTo_0));

idx = find(resize(dI_data_real,dateFrom:dateTo_0)<s.cases_min & resize(dI_data_reported,dateFrom:dateTo_0)<s.cases_min); %#ok<MXFND> % & delta<1-s.ratio_threshold); %#ok<MXFND>
idx = dateFrom:max(idx);
dI_data_real(idx) = dI_data_reported(idx);             
X_o(idx) = dI_data_reported_old(idx);       X_o(dateFrom:min(idx)) = dI_data_reported_old(dateFrom:min(idx));
X_y(idx) = dI_data_reported_old(idx);       X_y(dateFrom:min(idx)) = dI_data_reported_young(dateFrom:min(idx));
delta = method_params(resize(dI_data_reported,dateFrom:dateTo_0)./resize(dI_data_real,dateFrom:dateTo_0));
X = dI_data_real;

obs_ratio_adj(dateFrom:dateTo_0) = smooth_series(delta*s.obs_ratio);

sa = struct;
sa.Xs = (1-sigma(1)).*X;
sa.Xo = X_o;
sa.Xy = X_y;
sa.Xa = X-sa.Xs;
sa.dIa_data_reported = dI_data_reported.*sigma;
sa.dIs_data_reported = dI_data_reported-sa.dIa_data_reported;
sa.loss_a = method_params(sa.Xa-sa.dIa_data_reported);
sa.loss_s = sa.Xs-sa.dIs_data_reported; idx = find(sa.loss_s<0); sa.loss_s(idx) = 0; %#ok<FNDSB>
sa.loss_s = method_params(sa.loss_s);
sa.loss_o = method_params(sa.Xo-dI_data_reported_old);
sa.loss_y = method_params(sa.Xy-dI_data_reported_young);

% store params
p.omega = omega;
p.p_T_death = p_T_death;
p.T_rec = T_rec;
p.x_rec = x_rec;
p.p_T_rec = p_T_rec;
p.lambda = lambda;
p.p_T_hosp = p_T_hosp;
p.p_T_sick = p_T_sick;
p.T_delay = T_delay;
p.T_test_to_result = T_test_to_result;
p.x_sick = x_sick;
p.T_sick = T_sick;
p.x_shift = xs;
p.p_T_shift = p_T_shift;
p.T_shift = T_shift;
p.rho = rho_real;
p.varsigma = varsigma;
p.omega_o = s.omega_o.*(gamma_hd);
p.omega_y = s.omega_y.*(gamma_hd);

    function [y,ys] = adjust_series(x) %#ok<DEFNU>
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

    function [x] = adjust_tail(x,k) %#ok<DEFNU>
        dx = x(T-k)-x(T-k-1);
        x(T-k+1) = x(T-k)+2/3*dx;
        x(T-k+2) = x(T-k+1)+1/3*dx;
        for j=3:k
            x(T-k+j) = x(T-k+j-1)+1/3*1/(j-1)*dx;
        end        
    end

    function [x] = extend_tail(x,k) %#ok<DEFNU>
        L = length(x);
        dx = x(L)-x(L-1);
        x(L+1) = x(L)+2/3*dx;
        x(L+2) = x(L+1)+1/3*dx;
        for j=3:k
            x(L+j) = x(L+j-1)+1/3*1/(j-1)*dx;
        end        
    end

    function [x] = get_wa(weight,Z,alpha,idxFrom)
        sz = size(weight);
        weight = weight(idxFrom+1:end,:);
        alpha = alpha(idxFrom+1:end,:);
        k = sz(2);k0 = sz(1);
        t = length(Z)-idxFrom;
        phi = alpha./(1:k);
        if k0==1
            W = repmat(weight(end:-1:1),t,1);
            A = repmat(phi(end:-1:1),t,1);
        else
            W = weight(:,end:-1:1);
            A = phi(:,end:-1:1);
        end
        J = repmat(1:k,t,1)+repmat((0:t-1)',1,k);
        L = 0*(k-1)+repmat((1:t)',1,k);
        Weight_mat = sparse(L,J,W);
        Alpha_mat = sparse(L,J,A);
        Weight_mat = Weight_mat./sum(Weight_mat,2);
        x = (Weight_mat.*Alpha_mat)*Z(end-t-k+2:end);
    end

    function [x] = get_wa_inv(weight,zvec,x0,alpha,idxFrom)
        sz = size(weight);
        weight = weight(idxFrom+1:end,:);
        alpha = alpha(idxFrom+1:end,:);
        k = sz(2);k0 = sz(1);
        t = length(zvec)-idxFrom;
        phi = alpha./(1:k);
        if k0==1
            W = repmat(weight(k:-1:1),t,1);
            A = repmat(phi(k:-1:1),t,1);
        else
            W = weight(:,(k:-1:1));
            A = phi(:,(k:-1:1));
        end        
        J = repmat(1:k,t,1)+repmat((0:t-1)',1,k);
        L = repmat((1:t)',1,k);
        Weight_mat = sparse(L(:),J(:),W(:));
        Alpha_mat = sparse(L(:),J(:),A(:));
        x = zeros(t+k-1,1);
        x(1:k-1) = x0(1:idxFrom-1);
        Weight_mat = Weight_mat./sum(Weight_mat,2);
        function [d] = solve_lineqn(xx0)
            d=(Weight_mat.*Alpha_mat)*xx0 - zvec(end-t+1:end);
        end
        x = fsolve(@solve_lineqn,x,optimoptions('fsolve','Display','off','Algorithm','Levenberg-Marquardt'));
    end

    function [y] = extend(x,t0)
        xlen = length(x);
        z = x(1)+zeros(xlen+t0,1);
        z(t0+1:end) = x;
        y = method_params(z);
    end

end