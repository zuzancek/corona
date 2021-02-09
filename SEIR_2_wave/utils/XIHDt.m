function [out] = XIHDt(x,p,s,init,dateFrom,dateTo,firstDate)

T = dateTo-dateFrom+1;
T_total = dateTo-firstDate+1;
shift_i = max(s.k_hosp,s.k_sick);
shift_h = max(s.k_death,s.k_rec);
shift = max(shift_i,shift_h);
burnin = dateFrom-firstDate;
I0 = init.I(firstDate-shift_i+2:dateTo);
H0 = init.H(firstDate-shift_h+2:dateTo);
D0 = init.D(firstDate:dateTo);

method_params = s.smoothing_method_params;
method_data = s.smoothing_method_data;

try
    rho = init.rho(end-T_total+1:end);
catch err %#ok<NASGU>
    rho = method_params(init.rho);
    rho = rho(dateFrom:dateTo);
end
try
    T_delay = init.T_delay(end-T_total+1:end);
catch err %#ok<NASGU>
    T_delay = zeros(T,1);
end 
rho = extend(rho,shift_i);
T_obs = T_delay+s.T_test.mean;
varsigma = method_params(init.varsigma);
varsigma = varsigma(dateFrom:dateTo);

% ********* initialization
X = method_data(x.NewCases(firstDate:dateTo));

% ********* arrays (key)
% new cases: Old/Young
X_o = X.*rho(end-T_total+1:end);
X_y = X-X_o;
% active cases: Young/old
I0 = method_data(I0);
I_o = zeros(T_total+shift,1); I_o(end-T_total-shift_i+1:end) = I0.*rho;
I_y = I0-I_o;
% hospitalizations: Old/Young
H0 = method_data(H0);
H_o = zeros(T_total+shift,1); H_o(end-T_total-shift_h+1:end) = H0.*p.ho_h(1);
H_y = H0-H_o;
% deaths: Old/Young
D0 = method_data(D0);
D_o = zeros(T_total+shift,1); D_o(1:shift) = D0(1)*varsigma(1);
D_y = D0-D_o;

% other arrays
d_I_H_o = zeros(T,1);  d_I_H_y = d_I_H_o; d_I_R_o = d_I_H_o; d_I_R_y = d_I_H_o;
d_H_D_o = d_I_H_o; d_H_D_y = d_I_H_o; d_H_R_o = d_I_H_o; d_H_R_y = d_I_H_o;

% ******* parameters
% hospital admission
k_hosp = p.k_hosp;      t_hosp = s.time_h;
pdf_ih_o = extend(p.pdf_ih_o,T_total+k_hosp-size(p.pdf_ih_o,1)); alpha_iho = s.eta_o.*pdf_ih_o./repmat(t_hosp,T_total+k_hosp,1);
pdf_ih_y = extend(p.pdf_ih_y,T_total+k_hosp-size(p.pdf_ih_y,1)); alpha_ihy = s.eta_y.*pdf_ih_y./repmat(t_hosp,T_total+k_hosp,1);
alpha_iho = alpha_iho(:,end:-1:1);
alpha_ihy = alpha_ihy(:,end:-1:1);
% recovery from sickness at home
k_sick = p.k_sick;      
[pdf_ir_o,time_ir] = create_weights(k_sick,T_total+k_sick-size(p.pdf_ir_o,1),'Gamma',(s.T_sick_o-T_obs).*T_sick_std^2,1./s.T_sick_std^2);
pdf_ir_y = create_weights(k_sick,T_total+k_sick-size(p.pdf_ir_y,1),'Gamma',(s.T_sick_y-T_obs).*T_sick_std^2,1./s.T_sick_std^2);
alpha_iro = (1-s.eta_o).*pdf_ir_o./repmat(time_ir,T_total+k_sick,1);    alpha_iro = alpha_iro(:,end:-1:1);
alpha_iry = (1-s.eta_y).*pdf_ir_y./repmat(time_ir,T_total+k_sick,1);    alpha_iry = alpha_iry(:,end:-1:1);
% death at hospital
k_death = p.k_death;      t_death = s.time_d;
pdf_hd_o = extend(p.pdf_hd_o,T_total+k_death-size(p.pdf_hd_o,1)); alpha_hdo = s.omega_o.*pdf_hd_o./repmat(t_death,T_total+k_death,1);
pdf_hd_y = extend(p.pdf_hd_y,T_total+k_death-size(p.pdf_hd_y,1)); alpha_hdy = s.omega_y.*pdf_hd_y./repmat(t_death,T_total+k_death,1);
alpha_hdo = alpha_hdo(:,end:-1:1);
alpha_hdy = alpha_hdy(:,end:-1:1);
% recovery at hospital
k_rec = p.k_rec;      t_rec = s.time_r;
pdf_hr_o = extend(p.pdf_hr_o,T_total+k_rec-size(p.pdf_hr_o,1)); alpha_hro = (1-s.omega_o).*pdf_hr_o./repmat(t_rec,T_total+k_rec,1);
pdf_hr_y = extend(p.pdf_hr_y,T_total+k_rec-size(p.pdf_hr_y,1)); alpha_hry = (1-s.omega_y).*pdf_hr_y./repmat(t_rec,T_total+k_rec,1);
alpha_hro = alpha_hro(:,end:-1:1);
alpha_hry = alpha_hry(:,end:-1:1);

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
% D(t+1) = D(t)+H_D(t);

% ******* Calculation
for t=burnin+1:T
    d_I_H_o(t) = dot(alpha_iho,I_o(t-k_hosp:t-1));     d_I_R_o(t) = dot(alpha_iro(t,:),I_o(t-k_sick:t-1));
    I_o(t) = I_o(t-1)+X_o(t-delay)-d_I_H_o(t)-d_I_R_o(t);
    d_I_H_y(t) = dot(alpha_ihy,I_y(t-k_hosp:t-1));     d_I_R_y(t) = dot(alpha_iry(t,:),I_y(t-k_sick:t-1));
    I_y(t) = I_y(t-1)+X_y(t-delay)-d_I_H_y(t)-d_I_R_y(t);
    d_H_D_o(t) = dot(alpha_hdo(t,:),H_o(t-k_death:t-1));    d_H_R_o(t) = dot(alpha_hro(t,:),H_o(t-k_rec:t-1));
    H_o(t) = H_o(t-1)+d_I_H_o(t)-d_H_D_o(t)-d_H_R_o(t);
    d_H_D_y(t) = dot(alpha_hdy(t,:),H_y(t-k_death:t-1));    d_H_R_y(t) = dot(alpha_hry(t,:),H_y(t-k_rec:t-1));
    H_y(t) = H_y(t-1)+d_I_H_y(t)-d_H_D_y(t)-d_H_R_y(t);
    D_o(t) = D_o(t-1)+d_H_D_o(t);
    D_y(t) = D_y(t-1)+d_H_D_y(t);
end

% store results
out = struct;
out.X_o = tseries(dateFrom:dateTo,X_o);
out.X_y = tseries(dateFrom:dateTo,X_y);
out.I_o = tseries(dateFrom:dateTo,I_o(burnin+1:end));
out.I_y = tseries(dateFrom:dateTo,I_y(burnin+1:end));
out.H_o = tseries(dateFrom:dateTo,H_o(burnin+1:end));
out.H_y = tseries(dateFrom:dateTo,H_y(burnin+1:end));
out.D_o = tseries(dateFrom:dateTo,D_o(burnin+1:end));
out.D_y = tseries(dateFrom:dateTo,D_y(burnin+1:end));
out.X = out.X_o+out.X_y;
out.I = out.I_o+out.I_y;
out.H = out.H_o+out.H_y;
out.D = out.D_o+out.D_y;

figure;
pp1=plot(smooth_series(out.X),'linewidth',1,'Color','g');hold on; plot(out.X,'Color',[0.5 0.5 0.5]);hold on;
% pp2=plot(smooth_series(out.I),'linewidth',1,'Color',[0.2 0.6 0.33]);hold on;plot(out.I,'Color',[0.5 0.5 0.5]);hold on;
pp3=plot(smooth_series(resize(init.H,dateFrom:dateTo)),'linewidth',2,'Color','r');hold on;plot(resize(init.H,dateFrom:dateTo),'Color',[0.5 0.5 0.5]);hold on;
pp4=plot(smooth_series(out.H),'linewidth',2,'Color','m');hold on;plot(out.H,'Color',[0.5 0.5 0.5]);hold on;
pp5=plot(smooth_series(resize(init.D,dateFrom:dateTo)),'linewidth',2,'Color','c');hold on;plot(resize(init.D,dateFrom:dateTo),'Color',[0.5 0.5 0.5]);hold on;
pp6=plot(smooth_series(out.D),'linewidth',2,'Color','k');hold on;plot(out.D,'Color',[0.5 0.5 0.5]);hold on;
grid on;
legend([pp1 pp3 pp4 pp5 pp6],{'X','H_rep', 'H', 'D_rep','D'});

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