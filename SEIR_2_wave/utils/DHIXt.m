function [X,I,obs_ratio_adj,sa,p] = DHIXt(x,data,s,dateFrom,dateTo,t0,~,params,delay)

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

r = set_yo_ratios();

% death
k_death = s.k_death;
omega_y = s.omega_y;         pdf_d_y = repmat(s.pdf_d_y',length(varsigma),1);
omega_o = s.omega_o;         pdf_d_o = repmat(s.pdf_d_o',length(varsigma),1);
% rho = repmat(varsigma./(1-varsigma),1,k_death);
pdf_d = r.rho_ho_h.*pdf_d_o+(1-r.rho_ho_h).*pdf_d_y;
omega = sum(pdf_d,2);        pdf_d = pdf_d./omega;                 time_d = s.time_d;
weight_d = pdf_d./repmat(time_d',length(varsigma),1);
% recovery
k_rec = s.k_rec;
pdf_r_y = repmat(s.pdf_r_y',length(varsigma),1);
pdf_r_o = repmat(s.pdf_r_o',length(varsigma),1);
pdf_r = r.rho_ro_r.*pdf_r_o+(1-r.rho_ro_r).*pdf_r_y;
omega_r = sum(pdf_r,2);        pdf_r = pdf_r./omega_r;             time_r = s.time_r;
weight_r = pdf_r./repmat(time_r',length(varsigma),1);
% hospital admission
k_hosp = s.k_hosp;
eta_y = s.eta_y;        eta_o = s.eta_o;
pdf_h_y = repmat(s.pdf_h_y',length(varsigma),1);
pdf_h_o = repmat(s.pdf_h_o',length(varsigma),1);
pdf_h = r.rho_io_i.*pdf_h_o+(1-r.rho_io_i).*pdf_h_y;
eta = sum(pdf_h,2);        pdf_h = pdf_h./eta;             time_h = s.time_h;
weight_h = pdf_h./repmat(time_h',length(varsigma),1);
% recovery from sickness
T_delay = extend(T_delay,length(theta)-length(T_delay));
k_sick = s.k_sick; time_s = s.time_s;
pdf_s_y = pdf('Gamma',repmat(time_s',length(varsigma),1),repmat((s.T_sick_y-T_delay)*s.T_sick_std^2,1,k_sick),1/s.T_sick_std^2+zeros(length(T_delay),k_sick)); 
pdf_s_y = pdf_s_y(2:s.k_sick+1)./sum(pdf_s_y(2:s.k_sick+1));
pdf_s_o = pdf('Gamma',repmat(time_s',length(varsigma),1),repmat((s.T_sick_o-T_delay)*s.T_sick_std^2,1,k_sick),1/s.T_sick_std^2+zeros(length(T_delay),k_sick));pdf_s_o = pdf_s_o(2:s.k_sick+1)./sum(pdf_s_o(2:s.k_sick+1));
pdf_s = r.rho_so_s.*pdf_s_o+(1-r.rho_so_s).*pdf_s_y;
pdf_s = pdf_s./sum(pdf_s,2);             
weight_s = pdf_s./repmat(time_s',length(varsigma),1);

%% time shift (death->illness)
% ks = s.t_shift_clin;        xs = 1:ks;
% pp(1) = 1./(mean(T_death_o).*mean(varsigma));pp(2) = 1./(mean(T_death_y).*mean(1-varsigma));
% % pp(3) = 1./(mean(T_hosp_o).*mean(theta));pp(4) = 1./(mean(T_hosp_y).*mean(1-theta));
% lmat = repmat(pp,length(pp),1)-pp'; lmat(lmat==0) = NaN;
% p_T_shift = prod(pp).*sum(exp(-pp'.*xs)./repmat(prod(lmat,2,'omitnan'),1,ks),1);
% T_shift_death = ceil(dot(p_T_shift,xs));      % mean shift
% [~,idx] = (max(theta));
% T_shift_hosp = ceil(dot(xs(1:k_hosp)',p_T_hosp(idx,:)));
% T_shift_sick = ceil(dot(xs(1:k_sick)',p_T_sick(idx,:)));
% T_shift = T_shift_hosp+0*T_shift_death;
% ******* Equations
% I(t) = I(t-1)+X(t)-I_H(t)-I_R(t);     
% H(t) = H(t-1)+I_H(t)-H_D(t)-H_R(t);
% D(t) = D(t-1)+H_D(t);

% initialization
dI_data = method_data(x.NewCases(dateFrom:dateTo-cut));
dI_data_all = method_data(x.NewCases(dateFrom:dateTo));
D = x.Deaths(firstData:dateTo)*data.D(dateFrom)/x.Deaths(dateFrom); 
D(tshift+1:end) = data.D(dateFrom:dateTo);
H = data.H(firstData:dateTo);
AC = method_data(x.ActiveCases(firstData-k_hosp+2:dateTo));

% calculation
HD = method_data(D(2:end)-D(1:end-1));  
hd = method_params(get_wa(pdf_d,H,omega,k_death));
gamma_hd =  method_params(extend(HD(k_death+1:end)./hd(1:end-1),k_death));
omega = extend(method_params(omega(1:end-1).*gamma_hd),1);
HR = method_data(extend(get_wa(pdf_r,H,1-omega,k_rec),k_rec));
IH = H(2:end)-H(1:end-1)+HR(1:end-1)+HD;
I = method_data(get_wa_inv(pdf_h(1:end-1,:),IH,AC,eta(1:end-1),k_hosp));
IR = method_data(extend(get_wa(pdf_s(1:end-2,:),I,eta(1:end-2),k_sick),k_sick));
X = method_data(I(2:end)-I(1:end-1)+IR+IH);
X = extend_tail(X,3);

% X_orig = X;
% kk = (1+adj/T_shift).^(0:T_shift-1); 
% X(end-T_shift+1:end) = X(end-T_shift+1:end).*kk';

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

    function [r] = set_yo_ratios()
        r.do_dy = varsigma./(1-varsigma);
        a_d_y = (s.omega_y/s.T_death_y_mean); a_d_o = (s.omega_o/s.T_death_o_mean);
        r.ho_hy = (a_d_y./a_d_o).*r.do_dy;
        r.ho_h = r.ho_hy./(1+r.ho_hy);
        r.rho_ho_h = repmat(r.ho_h,1,s.k_death); %
        a_r_y = ((1-s.omega_y)/s.T_rec_y_mean); a_r_o = ((1-s.omega_o)/s.T_rec_o_mean);
        r.ro_ry = (a_r_o./a_r_y).*r.ho_hy;
        r.ro_r = r.ro_ry./(1+r.ro_ry);
        r.rho_ro_r = repmat(r.ro_r,1,s.k_rec); %
        a_h_y = s.eta_y./s.T_hosp_y_mean; a_h_o = s.eta_o./s.T_hosp_o_mean;
        r.io_iy = a_h_y./a_h_o.*(a_d_o+a_r_o)./(a_d_y+a_r_y).*r.ho_hy;
        r.io_i = r.io_iy./(1+r.io_iy);
        r.rho_io_i = repmat(r.io_i,1,s.k_hosp); %
        a_s_y = (1-s.eta_y)./s.T_sick_y_mean; a_s_o = (1-s.eta_o)./s.T_sick_o_mean;
        r.so_sy = (a_s_o./a_s_y).*r.so_sy;
        r.so_s = r.so_sy./(1+r.so_sy);
        r.rho_so_s = repmat(r.so_s,1,s.k_sick); %
    end

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