function [X,res,p] = DHIXt(x,data,s,dateFrom,dateTo,t0,params,delay)

%% initialization
T = dateTo-dateFrom+1;
method_data = s.smoothing_method_data; 
method_params = s.smoothing_method_params;
tshift = s.firstData_offset;
firstData = -tshift+dateFrom;
cut = 0*params.cutoff;

varsigma = extend(double(resize(params.death_old_ratio,dateFrom:dateTo)),tshift);
rho = method_params(params.cases_old_ratio(firstData:dateTo));
sigma = method_params(params.asymp_ratio(dateFrom:dateTo));
psi = extend(double(resize(params.serious_cases_ratio,dateFrom:dateTo)),tshift);
scale = s.S_H_rate;

%% testing
% delay in testing (gradual)
T_delay = create_series(delay);
T_obs = T_delay+s.T_test.mean;

%% data processing
% initialization
dI_data = method_data(x.NewCases(dateFrom:dateTo-cut));
dI_data_all = method_data(x.NewCases(dateFrom:dateTo));
D = method_data(data.D(firstData:dateTo)); % cummulative datehs
H = method_data(data.H(firstData:dateTo)); % hospitalizations
R = method_data(data.R(firstData:end)); % discharges from hospital
r0 = set_yo_ratios_params();

%% params initialization
% death
k_death = s.k_death;
pdf_hd_y = repmat(s.pdf_d_y',length(varsigma),1);
pdf_hd_o = repmat(s.pdf_d_o',length(varsigma),1);
% pdf_hd = r0.rho_ho_h.*pdf_hd_o+(1-r0.rho_ho_h).*pdf_hd_y;
% pdf_hd = pdf_hd./sum(pdf_hd,2);             
% kappa_dd = (s.omega_y.*pdf_hd_y)./(s.omega_o.*pdf_hd_o);
omega_o = r0.omega_o;   omega_y = r0.omega_y;
% recovery
k_rec = s.k_rec;
% pdf_hr_y = repmat(s.pdf_r_y',length(varsigma),1);
pdf_hr_y_s = repmat(s.pdf_r_y_s',length(varsigma),1);
pdf_hr_y_m = repmat(s.pdf_r_y_m',length(varsigma),1);
obj_r_y_s = s.obj_r_y_s;
obj_r_y = s.obj_r_y;
% pdf_hr_y = psi.*pdf_hr_y_s+(1-psi).*pdf_hr_y_m;
T_rec_m = r0.rho_ho_h.*s.T_rec_o_m_mean+(1-r0.rho_ho_h).*s.T_rec_y_m_mean;
T_rec_s = r0.rho_ho_h.*s.T_rec_o_s_mean+(1-r0.rho_ho_h).*s.T_rec_y_s_mean;
% T1 = sum(repmat([0:k_rec],length(varsigma),1).*pdf_hr_y,2);
% pdf_hr_o = repmat(s.pdf_r_o',length(varsigma),1);
pdf_hr_o_s = repmat(s.pdf_r_o_s',length(varsigma),1);
pdf_hr_o_m = repmat(s.pdf_r_o_m',length(varsigma),1);
obj_r_o_s = s.obj_r_o_s;
obj_r_o = s.obj_r_o;
% pdf_hr_o = psi.*pdf_hr_o_s+(1-psi).*pdf_hr_o_m;
% T2 = sum(repmat([0:k_rec],length(varsigma),1).*pdf_hr_o,2);
% pdf_hr = r0.rho_hro_hr.*pdf_hr_o+(1-r0.rho_hro_hr).*pdf_hr_y;
% pdf_hr = pdf_hr./sum(pdf_hr,2);             
% hospital admission
k_hosp = s.k_hosp;
pdf_ih_y = repmat(s.pdf_h_y',length(varsigma),1);
pdf_ih_o = repmat(s.pdf_h_o',length(varsigma),1);
% pdf_ih = r0.rho_io_i.*pdf_ih_o+(1-r0.rho_io_i).*pdf_ih_y;
% pdf_ih = pdf_ih./sum(pdf_ih,2);             % time_ih = s.time_h;
eta_o = r0.eta_o;       eta_y = r0.eta_y;
% recovery from sickness, mild cases, no need of hospital care
k_sick = s.k_sick; 
% [pdf_ir_y,time_ir] = create_weights(k_sick,length(varsigma),'Gamma',(s.T_sick_y-T_obs)*s.T_sick_std^2,1./s.T_sick_std^2); 
% pdf_ir_o = create_weights(k_sick,length(varsigma),'Gamma',(s.T_sick_o-T_obs)*s.T_sick_std^2,1./s.T_sick_std^2);
% pdf_ir = r0.rho_iro_ir.*pdf_ir_o+(1-r0.rho_iro_ir).*pdf_ir_y;
% pdf_ir = pdf_ir./sum(pdf_ir,2);     

I_ini = method_data(x.ActiveCases(firstData-k_hosp+2:dateTo));
I_o_ini = extend(r0.io_i,length(I_ini)-length(r0.io_i)).*I_ini; 
I_y_ini = I_ini-I_o_ini;
H_ini = method_data(data.H(firstData-s.k_death+2:dateTo));
H_o_ini = extend(r0.rho_ho_h(:,1),length(H_ini)-length(r0.rho_ho_h(:,1))).*H_ini;       
H_y_ini = H_ini-H_o_ini;

% ******* Equations, both in O/Y version + aggregation
% I(t) = I(t-1)+X(t)-I_H(t)-I_R(t);     
% H(t) = H(t-1)+I_H(t)-H_D(t)-H_R(t);
% D(t) = D(t-1)+H_D(t);

% calculation
% hospital deaths (Y,O,agg)
HD = extend(D(2:end)-D(1:end-1),1); 
HD_o = HD.*varsigma;                
HD_y = HD-HD_o;
h_o = max(0,method_data(get_wa_inv(pdf_hd_o,HD_o,H_o_ini,omega_o,k_death+1)));
h_y = max(0,method_data(get_wa_inv(pdf_hd_y,HD_y,H_y_ini,omega_y,k_death+1)));
h = max(1,h_o+h_y);
H_o = method_data(h_o./h).*H; 
H_y = H-H_o;
kappa_d = method_params(h./H); 
omega_o = omega_o.*repmat(kappa_d,1,k_death+1);
omega_y = omega_y.*repmat(kappa_d,1,k_death+1);
% recovery at hospital
% H_o_ini(end-length(H_o)+1:end) = H_o;
% H_y_ini(end-length(H_y)+1:end) = H_y;
zeta_o = repmat(1-omega_o(:,1),1,k_rec+1);
zeta_y = repmat(1-omega_y(:,1),1,k_rec+1);
w = kappa_d.*scale;
pdf_hr_o = w.*pdf_hr_o_s+pdf_hr_o_m.*(1-w);
pdf_hr_y = w.*pdf_hr_y_s+pdf_hr_y_m.*(1-w);
HR_o = (extend(get_wa(pdf_hr_o,H_o,zeta_o,k_rec+1),k_rec));
HR_y = (extend(get_wa(pdf_hr_y,H_y,zeta_y,k_rec+1),k_rec));
HR = HR_o+HR_y;
kappa_r = method_params(R./HR);
HR_o = HR_o.*kappa_r;
HR_y = HR_y.*kappa_r;
HR = HR_o+HR_y;
% hospital admission
IH_o = (extend(H_o(2:end)-H_o(1:end-1)+HR_o(2:end)+HD_o(2:end),1));
IH_y = (extend(H_y(2:end)-H_y(1:end-1)+HR_y(2:end)+HD_y(2:end),1));
IH = IH_y+IH_o;
% active cases (true)
varrho = (kappa_d-1).*scale.*T_rec_s(:,1)./T_rec_m(:,1);
eta_o = eta_o.*(1-0*varrho);
eta_y = eta_y.*(1-0*varrho);
I_o = (get_wa_inv(pdf_ih_o,IH_o,I_o_ini,eta_o,k_hosp+1)); I_o = method_params(extend(I_o(1:end-1),1));
I_y = (get_wa_inv(pdf_ih_y,IH_y,I_y_ini,eta_y,k_hosp+1)); I_y = method_params(extend(I_y(1:end-1),1));
I = I_o+I_y;
% recovered at home (no hospital needed)
[pdf_ir_y,time_ir] = create_weights(k_sick,length(varsigma),'Gamma',(s.T_sick_y+0*varrho./kappa_r.*s.T_rec_y_m_mean-T_obs)*s.T_sick_std^2,1./s.T_sick_std^2); 
pdf_ir_o = create_weights(k_sick,length(varsigma),'Gamma',(s.T_sick_o+0*varrho./kappa_r.*s.T_rec_o_m_mean-T_obs)*s.T_sick_std^2,1./s.T_sick_std^2);
IR_o = (extend(get_wa(pdf_ir_o(:,:),I_o,1-eta_o,k_sick+1),k_sick));
IR_y = (extend(get_wa(pdf_ir_y(:,:),I_y,1-eta_y,k_sick+1),k_sick));
IR = IR_o+IR_y;
% inflow of new cases
X_o = method_data(I_o(2:end)-I_o(1:end-1)+IR_o(2:end)+IH_o(2:end));
X_y = method_data(I_y(2:end)-I_y(1:end-1)+IR_y(2:end)+IH_y(2:end));
X = X_o+X_y;

fcast_per = ceil(max(r0.T_hosp_mean));
Xts = smooth_series(X(tshift:end)); Xts = tseries(dateFrom:dateFrom+length(Xts)-1,Xts);
Xrts = (X(tshift:end)); Xrts = method_data(tseries(dateFrom:dateFrom+length(Xrts)-1,Xrts));
Orts = method_data(tseries(dateFrom:dateFrom+length(dI_data)-1,dI_data));
Ots = smooth_series(Orts);
Orts0 = tseries(dateFrom:dateFrom+length(dI_data_all)-1,dI_data_all);
Ots0 = smooth_series(Orts0);

figure;
pp1=bar(Orts,'FaceAlpha',0.85);hold on;
p=bar(resize(Xrts,dateFrom:dateTo-fcast_per-1),'FaceAlpha',0.7);
bar(resize(Xrts,dateTo-fcast_per:dateTo),'FaceAlpha',0.4,'FaceColor',p.FaceColor);
pp2=bar(mov_median(resize(params.h,dateFrom:dateTo)),'FaceAlpha',0.33,'FaceColor','k');
plot(Ots,'b','linewidth',2);
plot(resize(Xts,dateFrom:dateTo-fcast_per-1),'r','linewidth',2);
plot(resize(Xts,dateTo-fcast_per:dateTo),'r-.','linewidth',2);
plot(smooth_series(mov_median(resize(params.h,dateFrom:dateTo))),'k--','linewidth',1);
grid on;
legend([pp1 p pp2],{'New cases: officially reported','New cases: implied by hospitals', 'Patients at hospitals'});
title('New cases: reported vs. real');

Dis_rep = method_data(tseries(dateFrom:dateTo,data.R(end-length(Xts)+1:end)));
Dis = tseries(dateFrom:dateTo,HR(end-length(Xts)+1:end));
Adm_rep = method_data(tseries(dateFrom:dateTo,data.A(end-length(Xts)+1:end)));
Adm = tseries(dateFrom:dateTo,IH(end-length(Xts)+1:end));

figure;
subplot(2,1,1)
bar(Adm_rep);hold on;
plot(Adm,'linewidth',2); 
grid on;
title('Admissions');
legend({'Reported','Implied'});
subplot(2,1,2)
bar(Dis_rep);hold on;
plot(Dis,'linewidth',2); 
grid on;
title('Discharges');
legend({'Reported','Implied'});

figure;
subplot(2,1,1)
plot(tseries(dateTo-length(kappa_d)+1:dateTo,kappa_d),'linewidth',1);hold on;
plot(tseries(dateTo-length(kappa_r)+1:dateTo,kappa_r),'linewidth',1);hold on;
grid on;
legend({'exceeding deaths','exceeding recoveries'});
title('Parameters adjustment');
subplot(2,1,2)
plot(tseries(dateTo-length(w)+1:dateTo,w.*H),'linewidth',1);hold on;
plot(tseries(dateTo-length(psi)+1:dateTo,psi.*H),'linewidth',1);hold on;
grid on;
legend({'implied','observed'});
title('Serious cases');

%
Len = length(Xts)-fcast_per;
res = struct();
res.X_all = tseries(firstData:dateTo,[X;X(end)]);
res.X_smooth_all = smooth_series(res.X_all);
res.X_smooth = resize(Xts,dateFrom:dateFrom+Len);
res.X_smooth_total = resize(Xts,dateFrom:dateFrom+Len+fcast_per); 
res.X_forecast_smooth = resize(Xts,dateFrom+Len:dateFrom+Len);
res.X_raw = resize(Xrts,dateFrom:dateFrom+Len);
res.X_raw_total = resize(Xrts,dateFrom:dateFrom+Len+fcast_per);
res.X_forecast_raw = resize(Xrts,dateFrom+Len:dateFrom+Len);
res.X_rep_smooth = Ots;
res.X_rep_forecast_smooth = resize(Ots0,enddate(Ots)+1:dateTo);
res.X_rep_raw = Orts;

% adjust series endpoints and get ratio
% X = [X;X(end)]; X_o = [X_o;X_o(end)]; X_y = [X_y;X_y(end)];
rho_real = method_data(tseries(firstData:dateTo,[X_o;X_o(end)])./tseries(firstData:dateTo,[X;X(end)])); 
rho_real_smooth = method_params(rho_real);
X = X(tshift:end);
X_o = X_o(tshift:end);
X_y = X-X_o;
sigma = sigma(end-length(X)+1:end);
dateTo_X = dateFrom+length(X);
dateTo_R = dateFrom+length(dI_data);
dateTo_0 = min(dateTo_X,dateTo_R);
obs_ratio_adj = tseries(t0:dateTo_0,s.obs_ratio);
X = Xts; %tseries(dateFrom:dateTo_X,method_data(X));
dI_data_real = resize(X,firstData:dateTo_X);
dI_data_reported = Ots;
dI_data_reported_old = dI_data_reported.*rho(tshift:tshift+length(dI_data_reported)-1);
dI_data_reported_young = dI_data_reported-dI_data_reported_old;

% idx = find(resize(dI_data_real,dateFrom:dateTo_0)<s.cases_min & resize(dI_data_reported,dateFrom:dateTo_0)<s.cases_min); %#ok<MXFND> % & delta<1-s.ratio_threshold); %#ok<MXFND>
% idx = dateFrom:max(idx);
% dI_data_real(idx) = dI_data_reported(idx);             
% X_o(idx) = dI_data_reported_old(idx);       X_o(dateFrom:min(idx)) = dI_data_reported_old(dateFrom:min(idx));
% X_y(idx) = dI_data_reported_old(idx);       X_y(dateFrom:min(idx)) = dI_data_reported_young(dateFrom:min(idx));
delta = method_params(resize(dI_data_reported,dateFrom:dateTo_0)./resize(dI_data_real,dateFrom:dateTo_0));
% X = dI_data_real;

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

% results
res.I = I;
res.obs_ratio_adj = obs_ratio_adj;
res.sa = sa;
res.D = D; res.HD = HD; res.HD_o = HD_o; res.HD_y = HD_y; 
res.H = H; res.H_o = H_o; res.H_y = H_y;
res.HR = HR; res.HR_o = HR_o; res.HR_y = HR_y;
res.IH = IH; res.IH_o = IH_o; res.IH_y = IH_y;
res.IR = IR; res.IR_o = IR_o; res.IR_y = IR_y;
res.X = X; res.X_o = X_o; res.X_y = X_y;
% store params
ini.I = I_ini; ini.I_o = I_o_ini; ini.I_y = I_y_ini;
ini.H = H_ini; ini.H_o = H_o_ini; ini.H_y = H_y_ini;
p = struct();
p.ini = ini;
p.varsigma = varsigma;
p.varrho = varrho;
p.rho = rho_real;
p.rho_smooth = rho_real_smooth;
p.sigma = sigma;
p.burnin = tshift;
p.T_delay = T_delay;
p.T_obs = T_obs;
p.par = r0;
p.pdf_hd_y = pdf_hd_y;
p.pdf_hd_o = pdf_hd_o;
p.pdf_hr_y = pdf_hr_y;
p.pdf_hr_o = pdf_hr_o;
p.time_ir = time_ir;
p.pdf_ih_y = pdf_ih_y;
p.pdf_ih_o = pdf_ih_o;
p.pdf_ir_y = pdf_ir_y;
p.pdf_ir_o = pdf_ir_o;
p.omega_o = omega_o;
p.omega_y = omega_y;
p.zeta_o = zeta_o;
p.zeta_y = zeta_y;
p.eta_o = eta_o;
p.eta_y = eta_y;
p.kappa_r = kappa_r;

    function [r] = adjust_distributions(k) %#ok<DEFNU>
        pp = obj_r_y; 
        obj_r_y = pp;
        L = length(k);
        dd = makedist('Burr','alpha',pp.alpha.*k,'c',pp.c+zeros(L,1),'k',pp.k+zeros(L,1));
        r.pdf_r_y = pdf(dd,0:k_rec); r.pdf_r_y = r.pdf_r_y./sum(r.pdf_r_y,2);
        pp = obj_r_o; pp.B = pp.B.*k; obj_r_o = pp;
        r.pdf_r_o = pdf(obj_r_o,0:k_rec); r.pdf_r_o = r.pdf_r_o./sum(r.pdf_r_o,2);
        pp = obj_r_y_s; pp.omega = pp.omega.*k; obj_r_y_s = pp;
        r.pdf_r_y_s = pdf(obj_r_y_s,0:k_rec); r.pdf_r_y_s = r.pdf_r_y_s./sum(r.pdf_r_y_s,2);
        pp = obj_r_o_s; pp.omega = pp.omega.*k; obj_r_o_s = pp;
        r.pdf_r_o_s = pdf(obj_r_o_s,0:k_rec); r.pdf_r_o_s = r.pdf_r_o_s./sum(r.pdf_r_o_s,2);
    end

    function [r] = set_yo_ratios_params()
        r.hdo_hdy = varsigma./(1-varsigma);
        a_hd_y = (s.omega_y/s.T_death_y_mean); a_hd_o = (s.omega_o/s.T_death_o_mean);
        r.ho_hy = (a_hd_y./a_hd_o).*r.hdo_hdy;
        r.ho_h = r.ho_hy./(1+r.ho_hy);
        r.rho_ho_h = repmat(r.ho_h,1,s.k_death+1); 
        r.omega_o = s.omega_o+zeros(length(varsigma),s.k_death+1);
        r.omega_y = s.omega_y+zeros(length(varsigma),s.k_death+1);
        r.omega = s.omega_o.*r.rho_ho_h+s.omega_y.*(1-r.rho_ho_h);
        %
        a_hr_y = ((1-s.omega_y)/s.T_rec_y_mean); a_hr_o = ((1-s.omega_o)/s.T_rec_o_mean);
        r.hro_hry = (a_hr_o./a_hr_y).*r.ho_hy;
        r.hro_hr = r.hro_hry./(1+r.hro_hry);
        r.rho_hro_hr = repmat(r.hro_hr,1,s.k_rec+1); 
        %
        r.iho_ihy = (a_hd_o+a_hr_o)./(a_hd_y+a_hr_y).*r.ho_hy;
        a_ih_y = s.eta_y./s.T_hosp_y_mean; a_ih_o = s.eta_o./s.T_hosp_o_mean;
        r.io_iy = a_ih_y./a_ih_o.*r.iho_ihy;
        r.io_i = r.io_iy./(1+r.io_iy);
        r.rho_io_i = repmat(r.io_i,1,s.k_hosp+1); 
        r.eta = s.eta_o.*r.rho_io_i+(1-r.rho_io_i).*s.eta_y;
        r.eta_o = s.eta_o+zeros(length(varsigma),s.k_hosp+1);
        r.eta_y = s.eta_y+zeros(length(varsigma),s.k_hosp+1);
        r.T_hosp_mean = s.T_hosp_o_mean.*r.rho_io_i+s.T_hosp_y_mean.*(1-r.rho_io_i);
        r.T_hosp_mean = r.T_hosp_mean(:,1);
        %
        a_ir_y = (1-s.eta_y)./s.T_sick_y_mean; a_ir_o = (1-s.eta_o)./s.T_sick_o_mean;
        r.iro_iry = (a_ir_o./a_ir_y).*r.io_iy;
        r.iro_ir = r.iro_iry./(1+r.iro_iry);
        r.rho_iro_ir = repmat(r.iro_ir,1,s.k_sick+1); %
        r.xo_xy = (a_ih_o+a_ir_o)./(a_ih_y+a_ir_y).*r.io_iy;
        r.rho_real_xo_x = r.xo_xy./(1+r.xo_xy);
        r.T_sick_mean = s.T_sick_o_mean.*r.rho_io_i+s.T_sick_y_mean.*(1-r.rho_io_i)-T_obs;
        r.T_sick_mean = r.T_sick_mean(:,1);
    end

    function [pdf_x,pnt_x] = create_weights(pnts_num,T_num,type,mean_x,stdev_x)
        weights = 0:pnts_num;
        pdf_x = pdf(type,repmat(weights,T_num,1),repmat(mean_x,1,pnts_num+1),repmat(stdev_x,T_num,pnts_num+1));
        pdf_x = pdf_x./sum(pdf_x,2); 
        pnt_x = repmat(weights,T_num,1);
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

    function [x,x_mat] = get_wa(weight,Z,alpha,idxFrom)
        sz = size(weight);
        weight = weight(idxFrom:end,:);
        alpha = alpha(idxFrom:end,:);
        k = sz(2);k0 = sz(1);
        t = length(Z)-idxFrom+1;
        phi = alpha./repmat((0:k-1),t,1); phi(:,1)=0; 
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
        Z = Z(end-t-k+2:end);
        x = (Weight_mat.*Alpha_mat)*Z;
        zvec = repmat((1:t)',1,k)+repmat((0:k-1),t,1);
        Z_mat = Z(zvec);  
        W = W./sum(W,2);
        x_mat = (W.*A).*Z_mat;
    end

    function [x] = get_wa_inv(weight,zvec,x0,alpha,idxFrom)
        sz = size(weight);
        weight = weight(idxFrom:end,:);
        alpha = alpha(idxFrom:end,:);
        k = sz(2);k0 = sz(1);
        t = length(zvec)-idxFrom+1;        
        phi = alpha./repmat((0:k-1),t,1); phi(:,1)=0; 
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
        [xlen,xwid] = size(x);
        z = x(1,:)+zeros(xlen+t0,xwid);
        z(t0+1:end,:) = x;
        y = method_data(z);
    end

    function [y] = create_series(str)
        y = NaN+zeros(T,1); y(1) = 0;
        dlen = length(str.v);
        for ii=1:dlen
            y(str.at(ii)-dateFrom) = str.v(ii);
        end
        if dlen
            y(T) = str.v(end);
        end
        y = method_params(interp1(find(~isnan(y)),y(find(~isnan(y))),1:T)'); %#ok<FNDSB>
        y = extend(y,length(varsigma)-length(y));
    end

end