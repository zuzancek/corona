function [out_rep,out_full] = XIHDt(x,p,s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
burnin = s.firstData_offset;
firstData = -burnin+dateFrom;

T_total = dateTo-firstData+1;
shift_i = max(s.k_hosp,s.k_sick);
shift_h = max(s.k_death,s.k_rec);
shift = max(shift_i,shift_h);

I0 = init.I(firstData-shift_i+2:dateTo);
H0 = init.H(firstData-shift_h+2:dateTo);
D0 = init.D(firstData:dateTo);
IH0 = init.IH(end-T+1:end);
HR0 = init.HR(end-T+1:end);
omega_o = init.omega_o;
omega_y = init.omega_y;

method_params = s.smoothing_method_params;
method_data = s.smoothing_method_data;

try
    rho = init.rho(end-T_total+1:end);
catch err %#ok<NASGU>
    rho = method_params(init.rho);
    rho = rho(dateFrom:dateTo);
end
try
    varrho = init.varrho(end-T_total+1:end);
catch err %#ok<NASGU>
    varrho = zeros(T,1);
end
try
    T_delay = init.T_delay(end-T_total+1:end);
catch err %#ok<NASGU>
    T_delay = zeros(T,1);
end 
rho = extend(rho,shift_i);
varrho = extend(varrho,T_total-length(varrho));
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
% hospitalizations: Old/Young
H0 = extend(H0,T_total+shift-length(H0));
H_o = H0.*extend(p.par.ho_h,T_total+shift-length(p.par.ho_h));
H_y = H0-H_o;
% deaths: Old/Young
D0 = extend(D0,T_total+shift-length(D0));
D_o = zeros(T_total+shift,1); 
D_o(1:shift) = D0(1)*varsigma(1);
D_y = D0-D_o;

% other arrays
d_I_H_o = zeros(T_total,1);  d_I_H_y = d_I_H_o; d_I_R_o = d_I_H_o; d_I_R_y = d_I_H_o;
d_H_D_o = d_I_H_o; d_H_D_y = d_I_H_o; d_H_R_o = d_I_H_o; d_H_R_y = d_I_H_o;

% ******* parameters
% hospital admission
k_hosp = s.k_hosp;      t_hosp = s.time_h;
% pdf_ih_o = extend(p.pdf_ih_o,T_total+k_hosp-size(p.pdf_ih_o,1)); 
alpha_iho = s.eta_o.*(1-varrho).*p.pdf_ih_o./repmat(t_hosp,T_total,1);%repmat(t_hosp',T_total+k_hosp,1);
% pdf_ih_y = extend(p.pdf_ih_y,T_total+k_hosp-size(p.pdf_ih_y,1)); 
alpha_ihy = s.eta_y.*(1-varrho).*p.pdf_ih_y./repmat(t_hosp,T_total,1);
alpha_iho = alpha_iho(:,end:-1:1);
alpha_ihy = alpha_ihy(:,end:-1:1);
% recovery from sickness at home
k_sick = s.k_sick;      
[pdf_ir_o,time_ir] = create_weights(k_sick,T_total+0*k_sick,'Gamma',(s.T_sick_o+varrho.*s.T_rec_o_m_mean-T_obs).*s.T_sick_std^2,1./s.T_sick_std^2);
pdf_ir_y = create_weights(k_sick,T_total+0*k_sick,'Gamma',(s.T_sick_y+varrho.*s.T_rec_o_m_mean-T_obs).*s.T_sick_std^2,1./s.T_sick_std^2);
alpha_iro = (1-s.eta_o).*pdf_ir_o./time_ir;    alpha_iro = alpha_iro(:,end:-1:1);
alpha_iry = (1-s.eta_y).*pdf_ir_y./time_ir;    alpha_iry = alpha_iry(:,end:-1:1);
% death at hospital
k_death = s.k_death;      t_death = s.time_d;
% pdf_hd_o = extend(p.pdf_hd_o,T_total+k_death-size(p.pdf_hd_o,1)); 
alpha_hdo = omega_o.*p.pdf_hd_o./repmat(t_death,T_total+0*k_death,1);
% pdf_hd_y = extend(p.pdf_hd_y,T_total+k_death-size(p.pdf_hd_y,1)); 
alpha_hdy = omega_y.*p.pdf_hd_y./repmat(t_death,T_total+0*k_death,1);
alpha_hdo = alpha_hdo(:,end:-1:1);
alpha_hdy = alpha_hdy(:,end:-1:1);
% recovery at hospital
k_rec = s.k_rec;      t_rec = s.time_r;
% pdf_hr_o = extend(p.pdf_hr_o,T_total+k_rec-size(p.pdf_hr_o,1)); 
alpha_hro = (1-omega_o).*p.pdf_hr_o./repmat(t_rec,T_total+0*k_rec,1);
% pdf_hr_y = extend(p.pdf_hr_y,T_total+k_rec-size(p.pdf_hr_y,1)); 
alpha_hry = (1-omega_y).*p.pdf_hr_y./repmat(t_rec,T_total+0*k_rec,1);
alpha_hro = alpha_hro(:,end:-1:1);
alpha_hry = alpha_hry(:,end:-1:1);

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
% D(t+1) = D(t)+H_D(t);

% ******* Calculation
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

% store results
out_full = struct;
out_full.X_o = tseries(firstData:dateTo,X_o(end-T_total+1:end));
out_full.X_y = tseries(firstData:dateTo,X_y(end-T_total+1:end));
out_full.I_o = tseries(firstData:dateTo,I_o(end-T_total+1:end));
out_full.I_y = tseries(firstData:dateTo,I_y(end-T_total+1:end));
out_full.H_o = tseries(firstData:dateTo,H_o(end-T_total+1:end));
out_full.H_y = tseries(firstData:dateTo,H_y(end-T_total+1:end));
out_full.D_o = tseries(firstData:dateTo,D_o(end-T_total+1:end));
out_full.D_y = tseries(firstData:dateTo,D_y(end-T_total+1:end));
out_full.IH_y = tseries(firstData:dateTo,d_I_H_y);
out_full.IH_o = tseries(firstData:dateTo,d_I_H_o);
out_full.HR_y = tseries(firstData:dateTo,d_H_R_y);
out_full.HR_o = tseries(firstData:dateTo,d_H_R_o);
out_full.X = out_full.X_o+out_full.X_y;
out_full.I = out_full.I_o+out_full.I_y;
out_full.H = out_full.H_o+out_full.H_y;
out_full.D = out_full.D_o+out_full.D_y;
out_full.IH = out_full.IH_o+out_full.IH_y;
out_full.HR = out_full.HR_o+out_full.HR_y;

out_rep = out_full;
out_rep.X_o = resize(out_rep.X_o,dateFrom:dateTo);
out_rep.X_y = resize(out_rep.X_y,dateFrom:dateTo);
out_rep.I_o = resize(out_rep.I_o,dateFrom:dateTo);
out_rep.I_y = resize(out_rep.I_y,dateFrom:dateTo);
out_rep.H_o = resize(out_rep.H_o,dateFrom:dateTo);
out_rep.H_y = resize(out_rep.H_y,dateFrom:dateTo);
out_rep.D_o = resize(out_rep.D_o,dateFrom:dateTo);
out_rep.D_y = resize(out_rep.D_y,dateFrom:dateTo);
out_rep.IH_y = resize(out_rep.IH_y,dateFrom:dateTo);
out_rep.IH_o = resize(out_rep.IH_o,dateFrom:dateTo);
out_rep.HR_y = resize(out_rep.HR_y,dateFrom:dateTo);
out_rep.HR_o = resize(out_rep.HR_o,dateFrom:dateTo);
out_rep.X = out_rep.X_o+out_rep.X_y;
out_rep.I = out_rep.I_o+out_rep.I_y;
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