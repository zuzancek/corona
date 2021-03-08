function [X,res,p] = DHIXs(x,data,s,dateFrom,dateTo,t0,params,delay)

%% initialization
T = dateTo-dateFrom+1;
method_data = s.smoothing_method_data; 
method_params = s.smoothing_method_params;
tshift = s.firstData_offset;
firstData = -tshift+dateFrom;
cut = 0*params.cutoff;
try
    hcut = 0*params.hosp_cutoff;
catch err %#ok<NASGU>
    hcut = 4;
end

varsigma = extend(double(resize(params.death_old_ratio,dateFrom:dateTo)),tshift);
rho = method_params(params.cases_old_ratio(firstData:dateTo));
mu = method_params(params.hy_ratio(firstData:dateTo));
sigma = method_params(params.asymp_ratio(dateFrom:dateTo));
N = s.sim_num;

%% testing
% delay in testing (gradual)
T_delay = create_series(delay);
T_obs = T_delay+s.T_test.mean;

%% data processing
% initialization
dI_data = method_data(x.NewCases(dateFrom:dateTo-cut));
dI_data_all = method_data(x.NewCases(dateFrom:dateTo));
D = method_data(data.D_raw(firstData:dateTo)); % cummulative deaths
HD = (data.D_raw(firstData+1:dateTo)-data.D_raw(firstData:dateTo-1));
delta_HD = mean(HD(end-3*hcut:end-hcut));
HD = extend(adjust_tail(HD,hcut,delta_HD)',1);
H = (data.H(firstData:dateTo)); % hospitalizations
S = (data.S(firstData:end)); % serious cases in hospitals
r0 = set_yo_ratios_params();

%% params initialization
% ******* 1./ hospital
% A./ death
k_death = s.k_death;
pdf_hd_y = repmat(s.pdf_hd_y',length(varsigma),1);
pdf_hd_o = repmat(s.pdf_hd_o',length(varsigma),1);
omega_o = r0.omega_o;   omega_y = r0.omega_y;
% B./ recovery
k_rec = s.k_rec;
pdf_hr_y = repmat(s.pdf_hr_y',length(varsigma),1);
pdf_hr_o = repmat(s.pdf_hr_o',length(varsigma),1);
T_rec_y_rv = random(s.obj_hr_y, N,1);
T_rec_o_rv = random(s.obj_hr_o, N,1);

% ******** 2./ home
% A./ admission to hospital
k_hosp = s.k_hosp;
pdf_ih_y = repmat(s.pdf_ih_y',length(varsigma),1);
pdf_ih_o = repmat(s.pdf_ih_o',length(varsigma),1);
eta_o = r0.eta_o;   eta_y = r0.eta_y;
% B./ recovery
k_sick = s.k_sick;
% mu_o = r0.mu_o;   mu_y = r0.mu_y;
T_rec_i_std = s.T_sick_std;
T_rec_i_y = s.T_sick_y;
T_rec_i_o = s.T_sick_o;
[pdf_ir_y,time_ir] = create_weights(k_sick,length(varsigma),'Gamma',(T_rec_i_y-T_obs)*T_rec_i_std^2,1./T_rec_i_std^2); 
pdf_ir_o = create_weights(k_sick,length(varsigma),'Gamma',(T_rec_i_o-T_obs)*T_rec_i_std^2,1./T_rec_i_std^2);
% ******* 3./ Serious cases (ICU,ECMO,...) - separate submodel
% A./ admission to ICU
k_ser = s.k_ser;
pdf_is_y = repmat(s.pdf_is_y',length(varsigma),1);
pdf_is_o = repmat(s.pdf_is_o',length(varsigma),1);
theta_o = r0.theta_o;   theta_y = r0.theta_y;
% B./ death
pdf_sd_y = repmat(s.pdf_sd_y',length(varsigma),1);
pdf_sd_o = repmat(s.pdf_sd_o',length(varsigma),1);
omega_o_s = r0.omega_o_s;   omega_y_s = r0.omega_y_s;
% C./ recovery
pdf_sr_y = repmat(s.pdf_sr_y',length(varsigma),1);
pdf_sr_o = repmat(s.pdf_sr_o',length(varsigma),1);

% initialization
I_ini = method_data(x.ActiveCases(firstData-k_hosp+2:dateTo));
I_o_ini = extend(r0.io_i,length(I_ini)-length(r0.io_i)).*I_ini; 
I_y_ini = I_ini-I_o_ini;
H_ini = method_data(data.H(firstData-s.k_death+2:dateTo));
H_o_ini = extend(r0.rho_ho_h(:,1),length(H_ini)-length(r0.rho_ho_h(:,1))).*H_ini;       
H_y_ini = H_ini-H_o_ini;
S_ini = method_data(data.S(firstData-s.k_death+2:dateTo));
S_o_ini = extend(r0.rho_ho_h(:,1),length(S_ini)-length(r0.rho_ho_h(:,1))).*S_ini;       
S_y_ini = S_ini-S_o_ini;
M_ini = H_ini-S_ini;
M_o_ini = H_o_ini-S_o_ini;
M_y_ini = H_y_ini-S_y_ini;

%% ******* Equations, both in O/Y version + aggregation
% I(t) = I(t-1)+X(t)-I_M(t)-I_R(t);     
% M(t) = M(t-1)+I_M(t)-M_S(t)-M_R(t);
% S(t) = S(t-1)+M_S(t)-S_D(t)-S_R(t);
% D(t) = D(t-1)+S_D(t);

%% calculation
% **** Hospital:
H_y = mu.*H;
H_o = H-H_y;
HD_o_imp = (extend(get_wa(pdf_hd_o,H_o,omega_o,k_death+1),k_death));
HD_y_imp = (extend(get_wa(pdf_hd_y,H_y,omega_y,k_death+1),k_death));
% kappa_d>1 <=> more people die than expected (based on hospitalization data)
% indication of more seriou cases
HD_o = HD.*varsigma;                
HD_y = HD-HD_o;
kappa_d_y = method_params(HD_y./HD_y_imp); omega_y = omega_y.*kappa_d_y;
kappa_d_o = method_params(HD_o./HD_o_imp); omega_o = omega_o.*kappa_d_o;
% recovery
zeta_o = repmat(1-omega_o(:,1),1,k_rec+1);
zeta_y = repmat(1-omega_y(:,1),1,k_rec+1);
HR_o = (extend(get_wa(pdf_hr_o,H_o,zeta_o,k_rec+1),k_rec));
HR_y = (extend(get_wa(pdf_hr_y,H_y,zeta_y,k_rec+1),k_rec));
HR = HR_o+HR_y;
xx = get_wa_rnd(T_rec_y_rv,pdf_hr_y,H_y,zeta_y,k_rec+1);
% admission to hospital
IH_o = (extend(H_o(2:end)-H_o(1:end-1)+HR_o(2:end)+HD_o(2:end),1));
IH_y = (extend(H_y(2:end)-H_y(1:end-1)+HR_y(2:end)+HD_y(2:end),1));
IH = IH_o+IH_y;
% ***** serious cases
% shares
s_o = theta_o(:,1)./eta_o(:,1).*H_o;
s_y = theta_y(:,1)./eta_y(:,1).*H_y;
s_tot = max(1,(s_o+s_y));
S_o = s_o./s_tot.*S;
S_y = S-S_o;
kappa_s = method_params(S./s_tot); 
% theta_o = kappa_s.*theta_o;
% theta_y = kappa_s.*theta_y;
% kappa_s> 1 <=> larger proportion of hospitalised patients are in more serious
% conditions than expected
% deaths
omega_o_s = min(1,omega_o_s.*repmat(kappa_d_o,1,k_death+1));
omega_y_s = min(1,omega_y_s.*repmat(kappa_d_y,1,k_death+1));
SD_o = (extend(get_wa(pdf_sd_o,S_o,omega_o_s,k_death+1),k_death));
SD_y = (extend(get_wa(pdf_sd_y,S_y,omega_y_s,k_death+1),k_death));
SD = SD_o+SD_y;
% recovery
zeta_o_s = repmat(1-omega_o_s(:,1),1,k_rec+1);
zeta_y_s = repmat(1-omega_y_s(:,1),1,k_rec+1);
SR_o = (extend(get_wa(pdf_sr_o,S_o,zeta_o_s,k_rec+1),k_rec));
SR_y = (extend(get_wa(pdf_sr_y,S_y,zeta_y_s,k_rec+1),k_rec));
SR = SR_o+SR_y;
% admission
IS_o = method_data(extend(S_o(2:end)-S_o(1:end-1)+SR_o(2:end)+SD_o(2:end),1));
IS_y = method_data(extend(S_y(2:end)-S_y(1:end-1)+SR_y(2:end)+SD_y(2:end),1));
IS = IS_o+IS_y;
% ***** home
% shares
i_o = (get_wa_inv(pdf_is_o,IS_o,I_o_ini,theta_o,k_ser+1)); i_o = method_params(extend(i_o(1:end-1),1));
i_y = (get_wa_inv(pdf_is_y,IS_y,I_y_ini,theta_y,k_ser+1)); i_y = method_params(extend(i_y(1:end-1),1));
% i = i_o+i_y;
I_o = (get_wa_inv(pdf_ih_o,IH_o,I_o_ini,eta_o,k_hosp+1)); I_o = method_params(extend(I_o(1:end-1),1));
I_y = (get_wa_inv(pdf_ih_y,IH_y,I_y_ini,eta_y,k_hosp+1)); I_y = method_params(extend(I_y(1:end-1),1));
% admission
kappa_h_o = method_params(max(1,I_o)./max(1,i_o)); eta_o = eta_o.*kappa_h_o;
kappa_h_y = method_params(max(1,I_y)./max(1,i_y)); eta_y = eta_y.*kappa_h_y;
I_o = (get_wa_inv(pdf_ih_o,IH_o,I_o_ini,eta_o,k_hosp+1)); I_o = method_params(extend(I_o(1:end-1),1));
I_y = (get_wa_inv(pdf_ih_y,IH_y,I_y_ini,eta_y,k_hosp+1)); I_y = method_params(extend(I_y(1:end-1),1));
I = I_o+I_y;
% recovery
mu_o = repmat(1-eta_o(:,1),1,k_sick+1);
mu_y = repmat(1-eta_y(:,1),1,k_sick+1);
IR_o = (extend(get_wa(pdf_ir_o(:,:),I_o,mu_o,k_sick+1),k_sick));
IR_y = (extend(get_wa(pdf_ir_y(:,:),I_y,mu_y,k_sick+1),k_sick));
IR = IR_o+IR_y;
% inflow of new cases
X_o = method_data(I_o(2:end)-I_o(1:end-1)+IR_o(2:end)+IH_o(2:end));
X_y = method_data(I_y(2:end)-I_y(1:end-1)+IR_y(2:end)+IH_y(2:end));
X = X_o+X_y;

%% aggregations
M = H-S; M_o = H_o-S_o; M_y = H_y-S_y;
MR = HR-SR; MR_o = HR_o-SR_o; MR_y = HR_y-SR_y;
MD = HD-SD; MD_o = HD_o-SD_o; MD_y = HD_y-SD_y;

fcast_per = ceil(max(r0.T_hosp_mean));
Xts = smooth_series(X(tshift:end)); Xts = tseries(dateFrom:dateFrom+length(Xts)-1,Xts);
Xrts = (X(tshift:end)); Xrts = method_data(tseries(dateFrom:dateFrom+length(Xrts)-1,Xrts));
Orts = method_data(tseries(dateFrom:dateFrom+length(dI_data)-1,dI_data));
Ots = smooth_series(Orts);
Orts0 = tseries(dateFrom:dateFrom+length(dI_data_all)-1,dI_data_all);
Ots0 = smooth_series(Orts0);

% figure;
% pp1=bar(Orts,'FaceAlpha',0.85);hold on;
% p=bar(resize(Xrts,dateFrom:dateTo-fcast_per-1),'FaceAlpha',0.7);
% bar(resize(Xrts,dateTo-fcast_per:dateTo),'FaceAlpha',0.4,'FaceColor',p.FaceColor);
% pp2=bar(mov_median(resize(params.h,dateFrom:dateTo)),'FaceAlpha',0.5,'FaceColor',[0.5 0.5 0.5]);
% pp3=bar(mov_median(resize(params.s,dateFrom:dateTo)),'FaceAlpha',0.5,'FaceColor','k');
% plot(Ots,'b','linewidth',2);
% plot(resize(Xts,dateFrom:dateTo-fcast_per-1),'r','linewidth',2);
% plot(resize(Xts,dateTo-fcast_per:dateTo),'r-.','linewidth',2);
% plot(smooth_series(mov_median(resize(params.h,dateFrom:dateTo))),'--','Color',[0.25 0.25 0.25],'linewidth',1);
% plot(smooth_series(mov_median(resize(params.s,dateFrom:dateTo))),'k-.','linewidth',1);
% grid on;
% legend([pp1 p pp2 pp3],{'New cases: officially reported','New cases: implied by hospitals', 'Hospitalizations', 'Intensive Care'});
% title('New cases: reported vs. real');

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
res.fcast_per = fcast_per;

% adjust series endpoints and get ratio
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
X = Xts; 
dI_data_real = resize(X,firstData:dateTo_X);
dI_data_reported = Ots;
dI_data_reported_old = dI_data_reported.*rho(tshift:tshift+length(dI_data_reported)-1);
dI_data_reported_young = dI_data_reported-dI_data_reported_old;
delta_raw = resize(dI_data_reported,dateFrom:dateTo_0)./resize(dI_data_real,dateFrom:dateTo_0);
delta = method_params(delta_raw);
obs_ratio_adj(dateFrom:dateTo_0) = delta*s.obs_ratio;
obs_ratio_adj_raw = delta_raw.*s.obs_ratio;
obs_ratio_ideal = 0*obs_ratio_adj_raw+s.obs_ratio;

sa = struct;
sa.Xs = (1-sigma).*X;
sa.Xo = X_o;
sa.Xy = X_y;
sa.Xa = X-sa.Xs;
res.sympt_share_ideal = s.symp_ratio_obs+0*X;
res.sympt_share_real = sa.Xs./X;

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
res.obs_ratio_adj_raw = obs_ratio_adj_raw;
res.obs_ratio_ideal = obs_ratio_ideal;
res.rho_real = rho_real;
res.rho_real_smooth = rho_real_smooth;
res.rho_obs = rho;
res.sa = sa;
res.D = D; res.SD = SD; res.SD_o = SD_o; res.SD_y = SD_y; 
res.S = S; res.S_o = S_o; res.S_y = S_y;
res.SR = SR; res.SR_o = SR_o; res.SR_y = SR_y;
res.M = M; res.M_o = M_o; res.M_y = M_y;
res.MR = MR; res.MR_o = MR_o; res.MR_y = MR_y;
res.MD = MD; res.MD_o = MD_o; res.MD_y = MD_y;
res.H = H; res.H_o = H_o; res.H_y = H_y;
res.IH = IH; res.IH_o = IH_o; res.IH_y = IH_y;
res.IS = IS; res.IS_o = IS_o; res.IS_y = IS_y;
res.IR = IR; res.IR_o = IR_o; res.IR_y = IR_y;
res.HR = HR; res.HR_o = HR_o; res.HR_y = HR_y;
res.X = X; res.X_o = X_o; res.X_y = X_y;
% store params
ini.I = I_ini; ini.I_o = I_o_ini; ini.I_y = I_y_ini;
ini.H = H_ini; ini.H_o = H_o_ini; ini.H_y = H_y_ini;
ini.S = S_ini; ini.S_o = S_o_ini; ini.S_y = S_y_ini;
ini.M = M_ini; ini.M_o = M_o_ini; ini.M_y = M_y_ini;

p = struct();
p.ini = ini;
p.varsigma = varsigma;
p.rho = rho_real;
p.rho_obs = rho;
p.rho_smooth = rho_real_smooth;
p.sigma = sigma;
p.burnin = tshift;
p.T_delay = T_delay;
p.T_obs = T_obs;
p.par = r0;
p.pdf_sd_y = pdf_sd_y;
p.pdf_sd_o = pdf_sd_o;
p.pdf_sr_y = pdf_sr_y;
p.pdf_sr_o = pdf_sr_o;
p.pdf_hr_y = pdf_hr_y;
p.pdf_hr_o = pdf_hr_o;
p.pdf_hd_y = pdf_hd_y;
p.pdf_hd_o = pdf_hd_o;
p.pdf_ir_y = pdf_ir_y;
p.pdf_ir_o = pdf_ir_o;
p.pdf_is_y = pdf_is_y;
p.pdf_is_o = pdf_is_o;
p.pdf_ih_y = pdf_ih_y;
p.pdf_ih_o = pdf_ih_o;
p.time_ir = time_ir;
p.omega_o = omega_o;
p.omega_y = omega_y;
p.omega_o_s = omega_o_s;
p.omega_y_s = omega_y_s;
p.zeta_o = zeta_o;
p.zeta_y = zeta_y;
p.eta_o = eta_o;
p.eta_y = eta_y;
p.mu_o = mu_o;
p.mu_y = mu_y;
p.mu = mu;
p.theta_y = theta_y;
p.theta_o = theta_o;
p.kappa_d_y = kappa_d_y;
p.kappa_d_o = kappa_d_o;
p.kappa_s = kappa_s;
p.kappa_h_o = kappa_h_o;
p.kappa_h_y = kappa_h_y;

    function [r] = set_yo_ratios_params()
        % death (h+s)
        r.omega_o = s.omega_o+zeros(length(varsigma),s.k_death+1);
        r.omega_y = s.omega_y+zeros(length(varsigma),s.k_death+1);
        r.omega_o_s = s.omega_o_s+zeros(length(varsigma),s.k_death+1);
        r.omega_y_s = s.omega_y_s+zeros(length(varsigma),s.k_death+1);
        % admission to ICU...
        r.theta_o = s.theta_o+zeros(length(varsigma),s.k_ser+1);
        r.theta_y = s.theta_y+zeros(length(varsigma),s.k_ser+1);
        r.vartheta_o = 1-s.theta_o+zeros(length(varsigma),s.k_rec+1);
        r.vartheta_y = 1-s.theta_y+zeros(length(varsigma),s.k_rec+1);
        % admission to hospital
        r.eta_o = s.eta_o+zeros(length(varsigma),s.k_hosp+1);
        r.eta_y = s.eta_y+zeros(length(varsigma),s.k_hosp+1);
        r.mu_o = 1-s.eta_o+zeros(length(varsigma),s.k_sick+1);
        r.mu_y = 1-s.eta_y+zeros(length(varsigma),s.k_sick+1); 
        r.hdo_hdy = varsigma./(1-varsigma);
        a_hd_y = (s.omega_y/s.T_death_y_mean); a_hd_o = (s.omega_o/s.T_death_o_mean);
        r.ho_hy = (a_hd_y./a_hd_o).*r.hdo_hdy;
        r.ho_h = r.ho_hy./(1+r.ho_hy);
        r.rho_ho_h = repmat(r.ho_h,1,s.k_death+1); 
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
        r.T_hosp_mean = s.T_hosp_o_mean.*r.rho_io_i+s.T_hosp_y_mean.*(1-r.rho_io_i);
        r.T_hosp_mean = r.T_hosp_mean(:,1);
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

    function [x] = adjust_tail(x,k,xt)
        x(end) = xt;
        x(end-k:end-1) = NaN;
        x = interp1(find(~isnan(x)),x(find(~isnan(x))),1:length(x)); %#ok<FNDSB>
    end

    function []=get_wa_rnd(tvec,tgrid,Z,alpha,idxFrom)
        tmax = tgrid(end);
        tidx = min(floor(tvec),tmax-1);
        tlam = tvec-tidx;
        tvecinv = 1./tvec;
        Z = Z(idxFrom:end);
        len = length(Z);
        zidx = -tidx+(1:len);
        zmat = tlam.*Z(zidx+1)+(1-tlam).*Z(zidx);
        z = alpha*tvecinv.*zmat;
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

%     function [x] = get_rv(y)
%         shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
%         L = length(shape0);
%         shape0_vec = repmat(shape0,N,1);
%         scale0_vec = scale0*ones(N,L);
%         x = gamrnd(shape0_vec,scale0_vec);
%     end

end