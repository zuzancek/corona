function [X,I,obs_ratio_adj,sa,p] = DHIXt(x,data,s,dateFrom,dateTo,t0,~,params,delay)

%% initialization
T = dateTo-dateFrom+1;
method_data = s.smoothing_method_data; 
method_params = s.smoothing_method_params;
tshift = s.firstData_offset;
firstData = -tshift+dateFrom;
cut = 0*params.cutoff;

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

r0 = set_yo_ratios_params();

% death
k_death = s.k_death;
pdf_hd_y = repmat(s.pdf_d_y',length(varsigma),1);
pdf_hd_o = repmat(s.pdf_d_o',length(varsigma),1);
pdf_hd = r0.rho_ho_h.*pdf_hd_o+(1-r0.rho_ho_h).*pdf_hd_y;
pdf_hd = pdf_hd./sum(pdf_hd,2);             omega = r0.omega;         % time_hd = s.time_d;
gamma_hd = (s.omega_y.*pdf_hd_y)/(s.omega_o.*pdf_hd_o);
% weight_hd = pdf_hd./repmat(time_hd',length(varsigma),1);
% recovery
k_rec = s.k_rec;
pdf_hr_y = repmat(s.pdf_r_y',length(varsigma),1);
pdf_hr_o = repmat(s.pdf_r_o',length(varsigma),1);
pdf_hr = r0.rho_hro_hr.*pdf_hr_o+(1-r0.rho_hro_hr).*pdf_hr_y;
pdf_hr = pdf_hr./sum(pdf_hr,2);             % time_hr = s.time_r;
% weight_hr = pdf_hr./repmat(time_hr',length(varsigma),1);
% hospital admission
k_hosp = s.k_hosp;
pdf_ih_y = repmat(s.pdf_h_y',length(varsigma),1);
pdf_ih_o = repmat(s.pdf_h_o',length(varsigma),1);
pdf_ih = r0.rho_io_i.*pdf_ih_o+(1-r0.rho_io_i).*pdf_ih_y;
pdf_ih = pdf_ih./sum(pdf_ih,2);             % time_ih = s.time_h;
eta = r0.eta;
% weight_ih = pdf_ih./repmat(time_ih',length(varsigma),1);
% recovery from sickness, mild cases, no need of hospital care
T_delay = extend(T_delay,length(varsigma)-length(T_delay));
T_obs = T_delay+s.T_test.mean;
k_sick = s.k_sick; 
[pdf_ir_y,time_ir] = create_weights(k_sick,length(varsigma),'Gamma',(s.T_sick_y-T_obs)*s.T_sick_std^2,1./s.T_sick_std^2); %#ok<ASGLU>
pdf_ir_o = create_weights(k_sick,length(varsigma),'Gamma',(s.T_sick_o-T_obs)*s.T_sick_std^2,1./s.T_sick_std^2);
pdf_ir = r0.rho_iro_ir.*pdf_ir_o+(1-r0.rho_iro_ir).*pdf_ir_y;
pdf_ir = pdf_ir./sum(pdf_ir,2);             
% weight_ir = pdf_ir./time_ir;

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
HD = extend(method_data(D(2:end)-D(1:end-1)),1);  
hd = method_params(get_wa(pdf_hd,H,omega,k_death));
gamma_hd =  method_params(extend(HD(k_death+1:end)./hd,k_death));
omega = repmat((method_params(omega(:,1).*gamma_hd)),1,k_death;
HR = method_data(extend(get_wa(pdf_hr,H,1-omega,k_rec),k_rec));
IH = extend(H(2:end)-H(1:end-1)+HR(2:end)+HD(2:end),1);
I = extend(method_data(get_wa_inv(pdf_ih,IH,AC,eta,k_hosp)),1);
IR = method_data(extend(get_wa(pdf_ir(:,:),I,1-eta,k_sick),k_sick));
X = method_data(I(2:end)-I(1:end-1)+IR(2:end)+IH(2:end));

Xts = smooth_series(X(tshift:end)); Xts = tseries(dateFrom:dateFrom+length(Xts)-1,Xts);
Xrts = (X(tshift:end)); Xrts = tseries(dateFrom:dateFrom+length(Xrts)-1,Xrts);
Orts = method_data(tseries(dateFrom:dateFrom+length(dI_data)-1,dI_data));
Ots = smooth_series(Orts);
Orts0 = tseries(dateFrom:dateFrom+length(dI_data_all)-1,dI_data_all);
Ots0 = smooth_series(Orts0);

figure;
bar(Orts);hold on;
bar(Xrts);
bar(mov_median(params.h));

Len = length(Xts);
p = struct();
p.X_smooth = resize(Xts,dateFrom:dateFrom+Len);
p.X_forecast_smooth = resize(Xts,dateFrom+Len:dateFrom+Len-1);
p.X_raw = resize(Xrts,dateFrom:dateFrom+Len);
p.X_forecast_raw = resize(Xrts,dateFrom+Len:dateFrom+Len-1);
p.X_rep_smooth = Ots;
p.X_rep_forecast_smooth = resize(Ots0,enddate(Ots)+1:dateTo);
p.X_rep_raw = Orts;
rho_real = method_params([r0.rho_real_xo_x(1);r0.rho_real_xo_x]);

% adjust series endpoints and get ratio
X = X(tshift:end);
rho_real = rho_real(end-length(X)+1:end);
sigma = sigma(end-length(X)+1:end);
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
% p.omega = omega;
% p.p_T_death = p_T_death;
% p.T_rec = T_rec;
% p.x_rec = x_rec;
% p.p_T_rec = p_T_rec;
% p.lambda = lambda;
% p.p_T_hosp = p_T_hosp;
% p.p_T_sick = p_T_sick;
% p.T_delay = T_delay;
% p.T_test_to_result = T_test_to_result;
% p.x_sick = x_sick;
% p.T_sick = T_sick;
% p.x_shift = xs;
% p.p_T_shift = p_T_shift;
% p.T_shift = T_shift;
% p.rho = rho_real;
% p.varsigma = varsigma;
% p.omega_o = s.omega_o.*(gamma_hd);
% p.omega_y = s.omega_y.*(gamma_hd);

    function [r] = set_yo_ratios_params()
        r.hdo_hdy = varsigma./(1-varsigma);
        a_hd_y = (s.omega_y/s.T_death_y_mean); a_hd_o = (s.omega_o/s.T_death_o_mean);
        r.ho_hy = (a_hd_y./a_hd_o).*r.hdo_hdy;
        r.ho_h = r.ho_hy./(1+r.ho_hy);
        r.rho_ho_h = repmat(r.ho_h,1,s.k_death); 
        r.omega = s.omega_o.*r.rho_ho_h+s.omega_y.*(1-r.rho_ho_h);
        %
        a_hr_y = ((1-s.omega_y)/s.T_rec_y_mean); a_hr_o = ((1-s.omega_o)/s.T_rec_o_mean);
        r.hro_hry = (a_hr_o./a_hr_y).*r.ho_hy;
        r.hro_hr = r.hro_hry./(1+r.hro_hry);
        r.rho_hro_hr = repmat(r.hro_hr,1,s.k_rec); 
        %
        r.iho_ihy = (a_hd_o+a_hr_o)./(a_hd_y+a_hr_y).*r.ho_hy;
        a_ih_y = s.eta_y./s.T_hosp_y_mean; a_ih_o = s.eta_o./s.T_hosp_o_mean;
        r.io_iy = a_ih_y./a_ih_o.*r.iho_ihy;
        r.io_i = r.io_iy./(1+r.io_iy);
        r.rho_io_i = repmat(r.io_i,1,s.k_hosp); 
        r.eta = s.eta_o.*r.rho_io_i+(1-r.rho_io_i).*s.eta_y;
        %
        a_ir_y = (1-s.eta_y)./s.T_sick_y_mean; a_ir_o = (1-s.eta_o)./s.T_sick_o_mean;
        r.iro_iry = (a_ir_o./a_ir_y).*r.io_iy;
        r.iro_ir = r.iro_iry./(1+r.iro_iry);
        r.rho_iro_ir = repmat(r.iro_ir,1,s.k_sick); %
        r.xo_xy = (a_ih_o+a_ir_o)./(a_ih_y+a_ir_y).*r.io_iy;
        r.rho_real_xo_x = r.xo_xy./(1+r.xo_xy);
    end

    function [pdf_x,pnt_x] = create_weights(pnts_num,T_num,type,mean_x,stdev_x)
        weights = 0:pnts_num;
        pdf_x = pdf(type,repmat(weights,T_num,1),repmat(mean_x,1,pnts_num+1),repmat(stdev_x,T_num,pnts_num+1));
        pdf_x = pdf_x./sum(pdf_x,2); 
        pdf_x = pdf_x(:,2:end);
        pnt_x = repmat(weights(2:end),T_num,1);
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