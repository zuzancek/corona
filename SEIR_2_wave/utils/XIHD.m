function [out_rep,out_full] = XIHD(x,p,s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
burnin = s.firstData_offset;
firstData = -burnin+dateFrom;

T_total = dateTo-firstData+1;
shift_i = max(s.k_hosp,s.k_sick);
shift_h = max(s.k_death,s.k_rec);
shift = max(shift_i,shift_h);

I0 = init.I(firstData-shift_i+2:dateTo);
H0 = init.H(firstData-shift_h+2:dateTo);
S0 = init.S(firstData-shift_h+2:dateTo);
D0 = init.D(firstData:dateTo);
IH0 = init.IH(end-T+1:end);
HR0 = init.HR(end-T+1:end);

method_params = s.smoothing_method_params;
method_data = s.smoothing_method_data;

try
    rho = init.rho(end-T_total+1:end);
catch err %#ok<NASGU>
    rho = method_params(init.rho);
    rho = rho(dateFrom:dateTo);
end
try
    kappa_s = init.kappa_s(end-T_total+1:end);
catch err %#ok<NASGU>
    kappa_s = ones(T,1);
end
try
    kappa_m = init.kappa_m(end-T_total+1:end);
catch err %#ok<NASGU>
    kappa_m = ones(T,1);
end
try
    T_delay = init.T_delay(end-T_total+1:end);
catch err %#ok<NASGU>
    T_delay = zeros(T,1);
end 
rho = extend(rho,shift_i);
kappa_s = extend(kappa_s,T_total-length(kappa_s));
kappa_m = extend(kappa_m,T_total-length(kappa_m));
T_obs = extend(T_delay+s.T_test.mean, length(rho)-length(T_delay)-shift_i);
varsigma = method_params(init.varsigma);
varsigma = extend(varsigma(dateFrom:dateTo),burnin);

% ********* initialization
X = method_data(x.NewCases(firstData:dateTo));

% ********* arrays (key)
% new cases: Old/Young
X_o = X.*rho(end-T_total+1:end);
X_y = X-X_o;
% active cases: Young/old
I0 = extend(I0,T_total+shift_i-length(I0));
I_o = zeros(T_total+shift,1); 
I_o(end-T_total-shift_i+1:end) = (I0.*rho(end-length(I0)+1:end));
I_y = extend(I0,shift_i)-I_o;
% hospitalizations: Old/Young, moderate cases (M) vs intensive care (S)
M0 = extend(H0-S0,T_total+shift-length(H0));
M_o = M0.*extend(p.par.ho_h,T_total+shift-length(p.par.ho_h));
M_y = M0-M_o;
S0 = extend(S0,T_total+shift-length(S0));
S_o = S0.*extend(p.par.ho_h,T_total+shift-length(p.par.ho_h));
S_y = S0-S_o;
% deaths: Old/Young
D0 = extend(D0,T_total+shift-length(D0));
D_o = zeros(T_total+shift,1); 
D_o(1:shift) = D0(1)*varsigma(1);
D_y = D0-D_o;

% other arrays
d_I_M_o = zeros(T_total,1);  d_I_M_y = d_I_M_o; d_I_R_o = d_I_M_o; d_I_R_y = d_I_M_o;
d_M_S_o = d_I_M_o; d_M_S_y = d_I_M_o; d_M_R_o = d_I_M_o; d_M_R_y = d_I_M_o;
d_S_D_o = d_I_M_o; d_S_D_y = d_I_M_o; d_S_R_o = d_I_M_o; d_S_R_y = d_I_M_o;

% ******* parameters
% ******* 1./ mild/asymptomatic cases, no hospital needed
% A./ hospital admission
k_hosp = s.k_hosp;      t_hosp = s.time_h;
alpha_imo = p.eta_o.*p.pdf_im_o./repmat(t_hosp,T_total,1);
alpha_imy = p.eta_y.*p.pdf_im_y./repmat(t_hosp,T_total,1);
alpha_imo = alpha_imo(:,end:-1:1);
alpha_imy = alpha_imy(:,end:-1:1);
% B./ recovery from sickness at home
k_sick = s.k_sick;      
[pdf_ir_o,time_ir] = create_weights(k_sick,T_total+0*k_sick,'Gamma',(s.T_sick_o-T_obs).*s.T_sick_std^2,1./s.T_sick_std^2);
pdf_ir_y = create_weights(k_sick,T_total+0*k_sick,'Gamma',(s.T_sick_y-T_obs).*s.T_sick_std^2,1./s.T_sick_std^2);
alpha_iro = (1-p.mu_o).*pdf_ir_o./time_ir;    alpha_iro = alpha_iro(:,end:-1:1);
alpha_iry = (1-p.mu_y).*pdf_ir_y./time_ir;    alpha_iry = alpha_iry(:,end:-1:1);
% ******** 2./ moderate cases, at hospital, no intensive case needed
% A./ admission to ICU...
k_ser = s.k_ser;  t_ser = s.time_ser;
pdf_ms_y = repmat(s.pdf_s_y',length(varsigma),1);
pdf_ms_o = repmat(s.pdf_s_o',length(varsigma),1);
theta_o = s.theta_o.*repmat(kappa_m,1,k_ser+1);
theta_y = s.theta_y.*repmat(kappa_m,1,k_ser+1);
alpha_mso = theta_o.*pdf_ms_o./repmat(t_ser,T_total,1);
alpha_msy = theta_y.*pdf_ms_y./repmat(t_ser,T_total,1);
alpha_mso = alpha_mso(:,end:-1:1);
alpha_msy = alpha_msy(:,end:-1:1);
% B./ recovery (at hospital, no IC)
k_rec = s.k_rec; t_rec = s.time_rec;
pdf_mr_y = repmat(s.pdf_mr_y',length(varsigma),1);
pdf_mr_o = repmat(s.pdf_mr_o',length(varsigma),1);
vartheta_o = 1-theta_o(:,1).*repmat(kappa_m,1,k_rec+1);
vartheta_y = 1-theta_y(:,1).*repmat(kappa_m,1,k_rec+1);
alpha_mro = vartheta_o.*pdf_mr_o./t_rec;    alpha_mro = alpha_mro(:,end:-1:1);
alpha_mry = vartheta_y.*pdf_mr_y./t_rec;    alpha_mry = alpha_mry(:,end:-1:1);
% ******* 3./ serious cases (hospital+Intensive care needed)
% A./ deaths
k_death = s.k_death; t_death = s.time_d;
pdf_sd_y = repmat(s.pdf_d_y',length(varsigma),1);
pdf_sd_o = repmat(s.pdf_d_o',length(varsigma),1);
omega_o = s.omega_o.*repmat(kappa_s,1,k_death+1);   
omega_y = s.omega_y.*repmat(kappa_s,1,k_death+1);
alpha_sdo = omega_o.*pdf_sd_o./repmat(t_death,T_total,1);
alpha_sdy = omega_y.*pdf_sd_y./repmat(t_death,T_total,1);
alpha_sdo = alpha_sdo(:,end:-1:1);
alpha_sdy = alpha_sdy(:,end:-1:1);
% B./ recovery
k_rec = s.k_rec;    t_rec = s.time_r;
pdf_sr_y = repmat(s.pdf_sr_y',length(varsigma),1);
pdf_sr_o = repmat(s.pdf_sr_o',length(varsigma),1);
zeta_o = repmat(1-omega_o(:,1),1,k_rec+1);
zeta_y = repmat(1-omega_y(:,1),1,k_rec+1);
alpha_sro = zeta_o.*pdf_sr_o./repmat(t_rec,T_total,1);
alpha_sry = zeta_y.*pdf_sr_y./repmat(t_rec,T_total,1);
alpha_sro = alpha_sro(:,end:-1:1);
alpha_sry = alpha_sry(:,end:-1:1);

% ******* Equations, both in O/Y version + aggregation
% I(t) = I(t-1)+X(t)-I_M(t)-I_R(t);     
% M(t) = M(t-1)+I_M(t)-M_S(t)-M_R(t);
% S(t) = S(t-1)+M_S(t)-S_D(t)-S_R(t);
% D(t) = D(t-1)+S_D(t);

% ******* Calculation
for t=1:T_total
    tt = t+shift;
    %
    d_I_M_o(t) = dot(alpha_imo(t,1:end-1),I_o(tt-k_hosp:tt-1));     
    d_I_R_o(t) = dot(alpha_iro(t,1:end-1),I_o(tt-k_sick:tt-1));
    I_o(tt) = I_o(tt-1)+X_o(t)-d_I_M_o(t)-d_I_R_o(t);
    %
    d_I_M_y(t) = dot(alpha_imy(t,1:end-1),I_y(tt-k_hosp:tt-1));     
    d_I_R_y(t) = dot(alpha_iry(t,1:end-1),I_y(tt-k_sick:tt-1));
    I_y(tt) = I_y(tt-1)+X_y(t)-d_I_M_y(t)-d_I_R_y(t);
    %
    d_M_S_o(t) = dot(alpha_mso(t,1:end-1),M_o(tt-k_ser:tt-1));     
    d_M_R_o(t) = dot(alpha_mro(t,1:end-1),M_o(tt-k_rec:tt-1));
    M_o(tt) = M_o(tt-1)+d_I_M_o(t)-d_M_S_o(t)-d_M_R_o(t);
    %
    d_M_S_y(t) = dot(alpha_msy(t,1:end-1),M_y(tt-k_ser:tt-1));     
    d_M_R_y(t) = dot(alpha_mry(t,1:end-1),M_y(tt-k_rec:tt-1));
    M_y(tt) = M_y(tt-1)+d_I_M_y(t)-d_M_S_y(t)-d_M_R_y(t);
    %
    d_S_D_o(t) = dot(alpha_sdo(t,1:end-1),S_o(tt-k_death:tt-1));    
    d_S_R_o(t) = dot(alpha_sro(t,1:end-1),S_o(tt-k_rec:tt-1));
    S_o(tt) = S_o(tt-1)+d_M_S_o(t)-d_S_D_o(t)-d_S_R_o(t);
    %
    d_S_D_y(t) = dot(alpha_sdy(t,1:end-1),S_y(tt-k_death:tt-1));    
    d_S_R_y(t) = dot(alpha_sry(t,1:end-1),S_y(tt-k_rec:tt-1));
    S_y(tt) = S_y(tt-1)+d_M_S_y(t)-d_S_D_y(t)-d_S_R_y(t);
    %
    D_o(tt) = D_o(tt-1)+d_S_D_o(t);
    D_y(tt) = D_y(tt-1)+d_S_D_y(t);
end

% store results
out_full = struct;
out_full.X_o = tseries(firstData:dateTo,X_o(end-T_total+1:end));
out_full.X_y = tseries(firstData:dateTo,X_y(end-T_total+1:end));
out_full.I_o = tseries(firstData:dateTo,I_o(end-T_total+1:end));
out_full.I_y = tseries(firstData:dateTo,I_y(end-T_total+1:end));
out_full.M_o = tseries(firstData:dateTo,M_o(end-T_total+1:end));
out_full.M_y = tseries(firstData:dateTo,M_y(end-T_total+1:end));
out_full.S_o = tseries(firstData:dateTo,S_o(end-T_total+1:end));
out_full.S_y = tseries(firstData:dateTo,S_y(end-T_total+1:end));
out_full.D_o = tseries(firstData:dateTo,D_o(end-T_total+1:end));
out_full.D_y = tseries(firstData:dateTo,D_y(end-T_total+1:end));
out_full.IH_y = tseries(firstData:dateTo,d_I_M_y);
out_full.IH_o = tseries(firstData:dateTo,d_I_M_o);
out_full.HR_y = tseries(firstData:dateTo,d_S_R_y+d_M_R_y);
out_full.HR_o = tseries(firstData:dateTo,d_S_R_o+d_M_R_o);
out_full.X = out_full.X_o+out_full.X_y;
out_full.I = out_full.I_o+out_full.I_y;
out_full.M = out_full.M_o+out_full.M_y;
out_full.S = out_full.S_o+out_full.S_y;
out_full.H_y = out_full.S_y+out_full.M_y;
out_full.H_o = out_full.S_o+out_full.M_o;
out_full.H = out_full.H_o+out_full.H_y;
out_full.D = out_full.D_o+out_full.D_y;
out_full.IH = out_full.IH_o+out_full.IH_y;
out_full.HR = out_full.HR_o+out_full.HR_y;

out_rep = out_full;
out_rep.X_o = resize(out_rep.X_o,dateFrom:dateTo);
out_rep.X_y = resize(out_rep.X_y,dateFrom:dateTo);
out_rep.I_o = resize(out_rep.I_o,dateFrom:dateTo);
out_rep.I_y = resize(out_rep.I_y,dateFrom:dateTo);
out_rep.M_o = resize(out_rep.M_o,dateFrom:dateTo);
out_rep.M_y = resize(out_rep.M_y,dateFrom:dateTo);
out_rep.S_o = resize(out_rep.S_o,dateFrom:dateTo);
out_rep.S_y = resize(out_rep.S_y,dateFrom:dateTo);
out_rep.D_o = resize(out_rep.D_o,dateFrom:dateTo);
out_rep.D_y = resize(out_rep.D_y,dateFrom:dateTo);
out_rep.IH_y = resize(out_rep.IH_y,dateFrom:dateTo);
out_rep.IH_o = resize(out_rep.IH_o,dateFrom:dateTo);
out_rep.HR_y = resize(out_rep.HR_y,dateFrom:dateTo);
out_rep.HR_o = resize(out_rep.HR_o,dateFrom:dateTo);
out_rep.X = out_rep.X_o+out_rep.X_y;
out_rep.M = out_rep.M_o+out_rep.M_y;
out_rep.S = out_rep.S_o+out_rep.S_y;
out_rep.H_y = out_rep.S_y+out_rep.M_y;
out_rep.H_o = out_rep.S_o+out_rep.M_o;
out_rep.H = out_rep.H_o+out_rep.H_y;
out_rep.D = out_rep.D_o+out_rep.D_y;
out_rep.IH = out_rep.IH_o+out_rep.IH_y;

figure;
pp3=bar(out_rep.H); hold on;
pp5=bar(out_rep.D,'FaceAlpha',0.5);hold on;
pp1=plot(smooth_series(out_rep.X),'linewidth',2,'Color',[0.75 0.75 0]);hold on; 
plot(out_rep.X,'Color',[0.5 0.5 0.5]);hold on;
pp2=plot((resize(init.H,dateFrom:dateTo)),'b','linewidth',2);hold on;
pp4=plot((resize(init.D,dateFrom:dateTo)),'r','linewidth',2);hold on;
grid on;
legend([pp1 pp2 pp3 pp4 pp5],{'New cases','Hospitalizations - reported', ...
    'Hospitalizations - implied', 'Deaths - reported, cummulative','Deaths - implied, cummulative'});

figure;
subplot(2,1,1)
pp1=bar(resize(init.A,dateFrom:dateTo),'FaceAlpha',0.5);hold on;
pp2=plot(out_rep.IH,'linewidth',2);
pp3=plot(tseries(dateFrom:dateTo,IH0),'linewidth',2,'Color',[1 1 1]*0.5);
grid on;
legend([pp1 pp2 pp3],{'Reported','Implied top-down (2x)','Implied bottom-up'});
title('Hospital admissions');

subplot(2,1,2)
pp1=bar(resize(init.R,dateFrom:dateTo),'FaceAlpha',0.5);hold on;
pp2=plot(out_rep.HR,'linewidth',2);
pp3=plot(tseries(dateFrom:dateTo,HR0),'linewidth',2,'Color',[1 1 1]*0.5);
grid on;
legend([pp1 pp2 pp3],{'Reported','Implied top-down (2x)','Implied bottom-up'});
title('Hospital discharges');

    function [y] = extend(x,t0)
        [xlen,xwid] = size(x);
        z = x(1,:)+zeros(xlen+t0,xwid);
        z(t0+1:end,:) = x;
        y = method_data(z);
    end

    function [pdf_x,pnt_x] = create_weights(pnts_num,T_num,type,mean_x,stdev_x)
        weights = 0:pnts_num;
        pdf_x = pdf(type,repmat(weights,T_num,1),repmat(mean_x,1,pnts_num+1),repmat(stdev_x,T_num,pnts_num+1));
        pdf_x = pdf_x./sum(pdf_x,2); 
        pnt_x = repmat(weights,T_num,1);
    end

end