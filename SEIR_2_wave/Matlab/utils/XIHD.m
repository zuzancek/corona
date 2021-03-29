function [out_rep,out_full] = XIHD(x,p,s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
burnin = s.firstData_offset;
firstData = -burnin+dateFrom;

T_total = dateTo-firstData+1;
shift_i = max([s.k_hosp,s.k_sick, s.k_ser]);
shift_h = max(s.k_death,s.k_rec);
shift = max(shift_i,shift_h);

I0 = init.I(firstData-shift_i+2:dateTo);
H0 = init.H(firstData-shift_h+2:dateTo);
S0 = init.S(firstData-shift_h+2:dateTo);
D0 = init.D(firstData:dateTo);

method_params = s.smoothing_method_params;
method_data = s.smoothing_method_data;

try
    rho = init.rho(end-T_total+1:end);
catch err %#ok<NASGU>
    rho = method_params(init.rho);
    rho = rho(dateFrom:dateTo);
end
try
    nu = double(init.nu);
catch err %#ok<NASGU>
    nu = ones(T,1);
end
try
    kappa_s = init.kappa_s(end-T_total+1:end);
catch err %#ok<NASGU>
    kappa_s = ones(T,1);
end
try
    kappa_h_y = init.kappa_h_y(end-T_total+1:end);
catch err %#ok<NASGU>
    kappa_h_y = ones(T,1);
end
try
    kappa_h_o = init.kappa_h_o(end-T_total+1:end);
catch err %#ok<NASGU>
    kappa_h_o = ones(T,1);
end
try
    kappa_d_y = init.kappa_d_y(end-T_total+1:end);
catch err %#ok<NASGU>
    kappa_d_y = ones(T,1);
end
try
    kappa_d_o = init.kappa_d_o(end-T_total+1:end);
catch err %#ok<NASGU>
    kappa_d_o = ones(T,1);
end
try
    adm_cut = init.adm_cut(end-T_total+1:end);
catch err %#ok<NASGU>
    adm_cut = zeros(T,1);
end 

rho = extend(rho,shift);
nu = extend(nu,T_total-length(nu));
kappa_h_o = extend(kappa_h_o,T_total-length(kappa_h_o));
kappa_h_y = extend(kappa_h_y,T_total-length(kappa_h_y));
kappa_s = extend(kappa_s,T_total-length(kappa_s));
kappa_d_y = extend(kappa_d_y,T_total-length(kappa_d_y));
kappa_d_o = extend(kappa_d_o,T_total-length(kappa_d_o));
adm_cut = extend(adm_cut,T_total-length(adm_cut));
varsigma = method_params(init.varsigma);
varsigma = extend(varsigma(dateFrom:dateTo),burnin);

% ********* initialization
X = method_data(x.NewCases(firstData:dateTo));

% ********* arrays (key)
% new cases: Old/Young
X_o = X.*rho(end-T_total+1:end);
X_y = X-X_o;
% active cases: Young/old
I0 = extend(I0,T_total+shift-length(I0));
I_o = 0*I0;
I_o(end-T_total-shift+1:end) = (I0.*rho(end-length(I0)+1:end));
I_y = I0-I_o;
J_o = I_o; J_y = I_y;
% hospitalizations: Old/Young, moderate cases (M) vs intensive care (S)
H0 = extend(H0,T_total+shift-length(H0));
H_o = H0.*extend(p.par.ho_h,T_total+shift-length(p.par.ho_h));
H_y = H0-H_o;
S0 = extend(S0,T_total+shift-length(S0));
S_o = S0.*extend(p.par.ho_h,T_total+shift-length(p.par.ho_h));
S_y = S0-S_o;
% deaths: Old/Young
D0 = extend(D0,T_total+shift-length(D0));
D_o = zeros(T_total+shift,1); 
D_o(1:shift) = D0(1)*varsigma(1);
D_y = D0-D_o;
F_y = D_y; F_o = D_o;

% params adjustment
eta_y = s.eta_y.*(1-.25*adm_cut).*kappa_h_y;
eta_o = s.eta_o.*kappa_h_o;
omega_o = s.omega_o.*kappa_d_o./nu;   
omega_y = s.omega_y.*kappa_d_y./nu;
theta_o = s.theta_o.*kappa_s./kappa_h_o;
theta_y = s.theta_y.*kappa_s./kappa_h_y;
omega_o_s = s.omega_o_s.*kappa_d_y./nu;
omega_y_s = s.omega_y_s.*kappa_d_o./nu;

% other arrays
d_I_H_o = zeros(T_total,1);  d_I_H_y = d_I_H_o; d_I_R_o = d_I_H_o; d_I_R_y = d_I_H_o;
d_J_S_o = d_I_H_o; d_J_S_y = d_I_H_o; d_J_R_o = d_I_H_o; d_J_R_y = d_I_H_o;
d_S_R_o = d_I_H_o; d_S_R_y = d_I_H_o; d_S_D_o = d_I_H_o; d_S_D_y = d_I_H_o; 
d_H_R_o = d_I_H_o; d_H_R_y = d_I_H_o; d_H_D_o = d_I_H_o; d_H_D_y = d_I_H_o;

% ******* parameters
% ******* 1./ mild/asymptomatic cases, no hospital needed
% A./ hospital admission
k_hosp = s.k_hosp;   
t_hosp = s.time_h;
alpha_iho = eta_o.*p.pdf_ih_o./repmat(t_hosp,T_total,1);
alpha_ihy = eta_y.*p.pdf_ih_y./repmat(t_hosp,T_total,1);
alpha_iho = alpha_iho(:,end:-1:1);
alpha_ihy = alpha_ihy(:,end:-1:1);
% B./ recovery from sickness at home
k_sick = s.k_sick;  
T_rec_i_y = s.T_sick_y-s.T_test.mean;
T_rec_i_o = s.T_sick_o-s.T_test.mean;
T_sick_std = s.T_sick_std;
[pdf_ir_y,time_ir] = create_weights(k_sick,T_total+0*k_sick,'Gamma',T_rec_i_y*T_sick_std^2+zeros(length(varsigma),1),1/(T_sick_std)^2);
pdf_ir_o = create_weights(k_sick,length(varsigma),'Gamma',T_rec_i_o*T_sick_std^2+zeros(length(varsigma),1),1/(T_sick_std)^2);
alpha_iro = (1-eta_o).*pdf_ir_o./time_ir;    alpha_iro = alpha_iro(:,end:-1:1);
alpha_iry = (1-eta_y).*pdf_ir_y./time_ir;    alpha_iry = alpha_iry(:,end:-1:1);
% ******** 2./ hospital care (general)
% A./ death
k_death = s.k_death; t_death = s.time_d;
pdf_hd_y = repmat(s.epdf_hd_y',length(varsigma),1);
pdf_hd_o = repmat(s.epdf_hd_o',length(varsigma),1);
% omega_o = s.omega_o.*repmat(kappa_d_o./nu,1,k_death+1);   
% omega_y = s.omega_y.*repmat(kappa_d_y./nu,1,k_death+1);
alpha_hdo = omega_o.*pdf_hd_o./repmat(t_death,T_total,1);
alpha_hdy = omega_y.*pdf_hd_y./repmat(t_death,T_total,1);
alpha_hdo = alpha_hdo(:,end:-1:1);
alpha_hdy = alpha_hdy(:,end:-1:1);
% B./ recovery
k_rec = s.k_rec; 
pdf_hr_y = repmat(s.epdf_hr_y',length(varsigma),1);
pdf_hr_o = repmat(s.epdf_hr_o',length(varsigma),1);
alpha_hro = (1-omega_o).*pdf_hr_o./s.time_r;    alpha_hro = alpha_hro(:,end:-1:1);
alpha_hry = (1-omega_y).*pdf_hr_y./s.time_r;    alpha_hry = alpha_hry(:,end:-1:1);
% ******* 3./ serious cases (hospital+Intensive care needed)
% A./ admission (ECMO, ICU, Ventilation)
k_ser = s.k_ser; t_ser = 0:k_ser;
alpha_jso = theta_o.*p.pdf_is_o./repmat(t_ser,T_total,1);
alpha_jsy = theta_y.*p.pdf_is_y./repmat(t_ser,T_total,1);
alpha_jso = alpha_jso(:,end:-1:1);
alpha_jsy = alpha_jsy(:,end:-1:1);
alpha_jro = (1-theta_o).*p.pdf_ir_o./time_ir; 
alpha_jry = (1-theta_y).*p.pdf_ir_y./time_ir; 
alpha_jro = alpha_jro(:,end:-1:1);
alpha_jry = alpha_jry(:,end:-1:1);
% B./ deaths
pdf_sd_y = repmat(s.epdf_sd_y',length(varsigma),1);
pdf_sd_o = repmat(s.epdf_sd_o',length(varsigma),1);
alpha_sdo = omega_o_s.*pdf_sd_o./repmat(t_death,T_total,1);
alpha_sdy = omega_y_s.*pdf_sd_y./repmat(t_death,T_total,1);
alpha_sdo = alpha_sdo(:,end:-1:1);
alpha_sdy = alpha_sdy(:,end:-1:1);
% C./ recovery
pdf_sr_y = repmat(s.epdf_sr_y',length(varsigma),1);
pdf_sr_o = repmat(s.epdf_sr_o',length(varsigma),1);
alpha_sro = (1-omega_o_s).*pdf_sr_o./s.time_r;
alpha_sry = (1-omega_y_s).*pdf_sr_y./s.time_r;
alpha_sro = alpha_sro(:,end:-1:1);
alpha_sry = alpha_sry(:,end:-1:1);

% ******* Equations, both in O/Y version + aggregation
% I(t) = I(t-1)+X(t)-I_M(t)-I_R(t);     
% M(t) = M(t-1)+I_M(t)-M_S(t)-M_R(t);
% S(t) = S(t-1)+M_S(t)-S_D(t)-S_R(t);
% D(t) = D(t-1)+S_D(t);

% ******* Calculation I.
for t=1:T_total
    tt = t+shift;
    %
    d_I_H_o(t) = dot(alpha_iho(t,1:end-1),I_o(tt-k_hosp:tt-1));     
    d_I_R_o(t) = dot(alpha_iro(t,1:end-1),I_o(tt-k_sick:tt-1));
    I_o(tt) = I_o(tt-1)+X_o(t)-d_I_H_o(t)-d_I_R_o(t);
    %
    d_I_H_y(t) = dot(alpha_ihy(t,1:end-1),I_y(tt-k_hosp:tt-1));     
    d_I_R_y(t) = dot(alpha_iry(t,1:end-1),I_y(tt-k_sick:tt-1));
    I_y(tt) = I_y(tt-1)+X_y(t)-d_I_H_y(t)-d_I_R_y(t);
    %
    d_H_D_o(t) = dot(alpha_hdo(t,1:end-1),H_o(tt-k_death:tt-1));     
    d_H_R_o(t) = dot(alpha_hro(t,1:end-1),H_o(tt-k_rec:tt-1));
    H_o(tt) = H_o(tt-1)+d_I_H_o(t)-d_H_D_o(t)-d_H_R_o(t);
    %    
    d_H_D_y(t) = dot(alpha_hdy(t,1:end-1),H_y(tt-k_death:tt-1));     
    d_H_R_y(t) = dot(alpha_hry(t,1:end-1),H_y(tt-k_rec:tt-1));
    H_y(tt) = H_y(tt-1)+d_I_H_y(t)-d_H_D_y(t)-d_H_R_y(t);    
    %
    D_o(tt) = D_o(tt-1)+d_H_D_o(t);
    D_y(tt) = D_y(tt-1)+d_H_D_y(t);
end

% ******* Calculation II.
for t=1:T_total
    tt = t+shift;
    %
    d_J_S_o(t) = dot(alpha_jso(t,1:end-1),J_o(tt-k_ser:tt-1));     
    d_J_R_o(t) = dot(alpha_jro(t,1:end-1),J_o(tt-k_sick:tt-1));
    J_o(tt) = J_o(tt-1)+X_o(t)-d_J_S_o(t)-d_J_R_o(t);
    %
    d_J_S_y(t) = dot(alpha_jsy(t,1:end-1),J_y(tt-k_ser:tt-1));     
    d_J_R_y(t) = dot(alpha_jry(t,1:end-1),J_y(tt-k_sick:tt-1));
    J_y(tt) = J_y(tt-1)+X_y(t)-d_J_S_y(t)-d_J_R_y(t);
    %    
    d_S_R_o(t) = dot(alpha_sro(t,1:end-1),S_o(tt-k_rec:tt-1));   
    d_S_D_o(t) = dot(alpha_sdo(t,1:end-1),S_o(tt-k_death:tt-1));
    S_o(tt) = S_o(tt-1)+d_J_S_o(t)-d_S_R_o(t)-d_S_D_o(t);
    %
    d_S_R_y(t) = dot(alpha_sry(t,1:end-1),S_y(tt-k_rec:tt-1));   
    d_S_D_y(t) = dot(alpha_sdy(t,1:end-1),S_y(tt-k_death:tt-1));
    S_y(tt) = S_y(tt-1)+d_J_S_y(t)-d_S_R_y(t)-d_S_D_y(t);
    %
    F_o(tt) = F_o(tt-1)+d_S_D_o(t);
    F_y(tt) = F_y(tt-1)+d_S_D_y(t);
end

%% process results
% store results
out_full = struct;
out_full.X_o = tseries(firstData:dateTo,X_o(end-T_total+1:end));
out_full.X_y = tseries(firstData:dateTo,X_y(end-T_total+1:end));
out_full.I_o = tseries(firstData:dateTo,I_o(end-T_total+1:end));
out_full.I_y = tseries(firstData:dateTo,I_y(end-T_total+1:end));
out_full.H_o = tseries(firstData:dateTo,H_o(end-T_total+1:end));
out_full.H_y = tseries(firstData:dateTo,H_y(end-T_total+1:end));
out_full.S_o = tseries(firstData:dateTo,S_o(end-T_total+1:end));
out_full.S_y = tseries(firstData:dateTo,S_y(end-T_total+1:end));
out_full.D_o = tseries(firstData:dateTo,D_o(end-T_total+1:end));
out_full.D_y = tseries(firstData:dateTo,D_y(end-T_total+1:end));
out_full.F_o = tseries(firstData:dateTo,F_o(end-T_total+1:end));
out_full.F_y = tseries(firstData:dateTo,F_y(end-T_total+1:end));
out_full.IH_y = tseries(firstData:dateTo,d_I_H_y);
out_full.IH_o = tseries(firstData:dateTo,d_I_H_o);
out_full.HR_y = tseries(firstData:dateTo,d_H_R_y);
out_full.HR_o = tseries(firstData:dateTo,d_H_R_o);
out_full.HD_y = tseries(firstData:dateTo,d_H_D_y);
out_full.HD_o = tseries(firstData:dateTo,d_H_D_o);
out_full.M_o = out_full.H_o-out_full.S_o;
out_full.M_y = out_full.H_y-out_full.S_y;
out_full.X = out_full.X_o+out_full.X_y;
out_full.I = out_full.I_o+out_full.I_y;
out_full.H = out_full.H_o+out_full.H_y;
out_full.M = out_full.M_o+out_full.M_y;
out_full.S = out_full.S_o+out_full.S_y;
out_full.D = out_full.D_o+out_full.D_y;
out_full.F = out_full.F_o+out_full.F_y;
out_full.IH = out_full.IH_o+out_full.IH_y;
out_full.HR = out_full.HR_o+out_full.HR_y;
out_full.HD = out_full.HD_o+out_full.HD_y;

out_rep = out_full;
out_rep.X_o = resize(out_rep.X_o,dateFrom:dateTo);
out_rep.X_y = resize(out_rep.X_y,dateFrom:dateTo);
out_rep.I_o = resize(out_rep.I_o,dateFrom:dateTo);
out_rep.I_y = resize(out_rep.I_y,dateFrom:dateTo);
out_rep.H_o = resize(out_rep.H_o,dateFrom:dateTo);
out_rep.H_y = resize(out_rep.H_y,dateFrom:dateTo);
out_rep.S_o = resize(out_rep.S_o,dateFrom:dateTo);
out_rep.S_y = resize(out_rep.S_y,dateFrom:dateTo);
out_rep.D_o = resize(out_rep.D_o,dateFrom:dateTo);
out_rep.D_y = resize(out_rep.D_y,dateFrom:dateTo);
out_rep.F_o = resize(out_rep.F_o,dateFrom:dateTo);
out_rep.F_y = resize(out_rep.F_y,dateFrom:dateTo);
out_rep.M_o = resize(out_rep.M_o,dateFrom:dateTo);
out_rep.M_y = resize(out_rep.M_y,dateFrom:dateTo);
out_rep.IH_y = resize(out_rep.IH_y,dateFrom:dateTo);
out_rep.IH_o = resize(out_rep.IH_o,dateFrom:dateTo);
out_rep.HR_y = resize(out_rep.HR_y,dateFrom:dateTo);
out_rep.HR_o = resize(out_rep.HR_o,dateFrom:dateTo);
out_rep.HD_y = resize(out_rep.HD_y,dateFrom:dateTo);
out_rep.HD_o = resize(out_rep.HD_o,dateFrom:dateTo);
out_rep.X = out_rep.X_o+out_rep.X_y;
out_rep.I = out_rep.I_o+out_rep.I_y;
out_rep.H = out_rep.H_o+out_rep.H_y;
out_rep.M = out_rep.M_o+out_rep.M_y;
out_rep.S = out_rep.S_o+out_rep.S_y;
out_rep.D = out_rep.D_o+out_rep.D_y;
out_rep.F = out_rep.F_o+out_rep.F_y;
out_rep.IH = out_rep.IH_o+out_rep.IH_y;
out_rep.HR = out_rep.HR_o+out_rep.HR_y;
out_rep.HD = out_rep.HD_o+out_rep.HD_y;

out_rep.varsigma = out_rep.HD_o./out_rep.HD;
out_rep.mu = out_rep.H_y./out_rep.H;
out_rep.eta_o = s.eta_o.*kappa_h_o;
out_rep.eta_y = s.eta_y.*kappa_h_y;

figure;
pp3=bar(out_rep.H); hold on;
pp5=bar(out_rep.D,'FaceAlpha',0.5);hold on;
pp7=bar(out_rep.S,'FaceColor',0.25*[1 1 1]); hold on;
pp1=plot(smooth_series(out_rep.X),'linewidth',2,'Color',[0.75 0.75 0]);hold on; 
plot(out_rep.X,'Color',[0.5 0.5 0.5]);hold on;
pp2=plot((resize(init.H,dateFrom:dateTo)),'b','linewidth',2);hold on;
pp6=plot((resize(init.S,dateFrom:dateTo)),'k','linewidth',2);hold on;
pp4=plot((resize(init.D,dateFrom:dateTo)),'r','linewidth',2);hold on;
grid on;
legend([pp1 pp2 pp3 pp4 pp5 pp6 pp7],{'New cases','Hospitalizations - reported', ...
    'Hospitalizations - implied', 'Deaths - reported, cummulative','Deaths - implied, cummulative', ...
    'Serious cases - reported', 'Serious cases - implied'});

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