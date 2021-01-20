function [out] = XIHDt(x,p,s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
I0 = init.I;     
H0 = init.H; 
D0 = init.D;

method_params = s.smoothing_method_params;
try
    rho = p.rho;
catch err
    rho = method_params(init.rho); rho = rho(dateFrom:dateTo);
end
varsigma = method_params(init.varsigma); varsigma = varsigma(dateFrom:dateTo);
     
eta_y = s.eta_y;                   
eta_o = s.eta_o;           

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);     
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
% D(t+1) = D(t)+H_D(t);

% initialization
method_data = s.smoothing_method_data;
X = method_data(x.NewCases(dateFrom:dateTo));
X_o = X.*rho;   X_y = X-X_o;
D0 = method_data(D0(dateFrom:dateTo));
D_o = zeros(T,1); D_o(1) = D0(1)*varsigma(1); 
D_y = zeros(T,1); D_y(1) = D0(1)-D_o(1); 
H0 = method_data(H0(dateFrom:dateTo));
H_o = zeros(T,1); H_o(1) = H0(1)*eta_o*rho(1)/(eta_o*s.old_share+eta_y*(1-s.old_share));
H_y = zeros(T,1); H_y(1) = H0(1)-H_o(1);
I0 = method_data(I0(dateFrom:dateTo));
I_o = zeros(T,1); I_o(1) = I0(1)*s.old_share;
I_y = zeros(T,1); I_y(1) = I0(1)-I_o(1);
d_I_H_o = zeros(T,1);  d_I_H_y = d_I_H_o; d_I_R_o = d_I_H_o; d_I_R_y = d_I_H_o;
d_H_D_o = d_I_H_o; d_H_D_y = d_I_H_o; d_H_R_o = d_I_H_o; d_H_R_y = d_I_H_o;

startFrom = max(s.k_hosp,s.k_sick,s.k_death,s.k_rec)+1;
T_total = T+startFrom;
T_delay = p.T_delay(end-T_total+1:end);

% parameters
% hospital admission
k_hosp = s.k_hosp;      x_hosp = (1:k_hosp)'; 
T_hosp_o = s.T_hosp_o;  p_T_hosp_o = pdf(s.T_hosp_pdf_type,x_hosp,1./T_hosp_o+zeros(k_hosp,1));
T_hosp_y = s.T_hosp_y;  p_T_hosp_y = pdf(s.T_hosp_pdf_type,x_hosp,1./T_hosp_y+zeros(k_hosp,1));
eta_o = s.eta_o;        alpha_iho = eta_o./x_hosp.*p_T_hosp_o;
eta_y = s.eta_y;        alpha_ihy = eta_y./x_hosp.*p_T_hosp_y;
% recovery from sickness at home
k_sick = s.k_sick;      x_sick = (1:k_sick)';       x_sick_mat = repmat(x_sick,T_total,1); T_sick_std2 = s.T_sick_std^2;
T_sick_o = s.T_sick_o-T_delay-p.Test_to_result;     T_sick_o_shape = T_sick_o*T_sick_std2; T_sick_scale = 1/T_sick_std2;
p_T_sick_o = pdf(s.T_sick_pdf_type,x_sick_mat,repmat(T_sick_o_shape,1,k_sick),repmat(T_sick_scale,1,k_sick));
T_sick_y = s.T_sick_y-T_delay-p.Test_to_result;     T_sick_y_shape = T_sick_y*T_sick_std2; 
p_T_sick_y = pdf(s.T_sick_pdf_type,x_sick_mat,repmat(T_sick_y_shape,1,k_sick),repmat(T_sick_scale,1,k_sick));
alpha_iro = (1-eta_o)./x_sick_mat.*p_T_sick_o;
alpha_iry = (1-eta_y)./x_sick_mat.*p_T_sick_y;
% death at hospital
k_death = s.k_death;        x_death = (1:k_death)';
T_death_o = s.T_death_o;    p_T_death_o = pdf(s.T_death_pdf_type,x_death,1./T_death_o+zeros(k_death,1));
T_death_y = s.T_death_y;    p_T_death_y = pdf(s.T_death_pdf_type,x_death,1./T_death_y+zeros(k_hosp,1));
omega_o = s.omega_o;        alpha_hdo = omega_o./x_death.*p_T_death_o;
omega_y = s.omega_y;        alpha_hdy = omega_y./x_death.*p_T_death_y;
% recovery at hospital
k_rec = s.k_rec;            x_rec = (1:k_rec)';                     T_rec_std2 = s.T_rec_std^2; T_rec_scale = 1/T_rec_std2+zeros(k_sick,1);
T_rec_o = s.T_rec_o;        T_rec_o_shape = T_rec_o*T_sick_std2;    p_T_rec_o = pdf(s.T_rec_pdf_type,x_rec,T_rec_o_shape+zeros(k_sick,1),T_rec_scale);
T_rec_y = s.T_rec_y;        T_rec_y_shape = T_rec_y*T_sick_std2;    p_T_rec_y = pdf(s.T_rec_pdf_type,x_rec,T_rec_y_shape+zeros(k_sick,1),T_rec_scale);
alpha_hro = (1-omega_o)./x_rec.*p_T_rec_o;
alpha_hry = (1-omega_y)./x_rec.*p_T_rec_y;

% calculation
for t=startFrom:T_total-1 
    d_I_H_o(t) = dot(alpha_iho,I_o(t-k_hosp+1:t));     d_I_R_o(t) = dot(alpha_iro(t,:),I_o(t-k_sick+1:t));
    I_o(t+1) = I_o(t)+X_o(t)-d_I_H_o(t)-d_I_R_o(t);
    d_I_H_y(t) = dot(alpha_ihy,I_y(t-k_hosp+1:t));     d_I_R_y(t) = dot(alpha_iry(t,:),I_y(t-k_sick+1:t));
    I_y(t+1) = I_y(t)+X_y(t)-d_I_H_y(t)-d_I_R_y(t);
    d_H_D_o(t) = dot(alpha_hdo,H_o(t-k_death+1:t));    d_H_R_o(t) = dot(alpha_hro,H_o(t-k_rec+1:t));
    H_o(t+1) = H_o(t)+d_I_H_o(t)-d_H_D_o(t)-d_H_R_o(t);
    d_H_D_y(t) = dot(alpha_hdy,H_y(t-k_death+1:t));    d_H_R_y(t) = dot(alpha_hry,H_y(t-k_rec+1:t));
    H_y(t+1) = H_y(t)+d_I_H_y(t)-d_H_D_y(t)-d_H_R_y(t);
    D_o(t+1) = D_o(t)+d_H_D_o(t);
    D_y(t+1) = D_y(t)+d_H_D_y(t);
end

% store results
out = struct;
out.X_o = tseries(dateFrom:dateTo,X_o);
out.X_y = tseries(dateFrom:dateTo,X_y);
out.I_o = tseries(dateFrom:dateTo,I_o);
out.I_y = tseries(dateFrom:dateTo,I_y);
out.H_o = tseries(dateFrom:dateTo,H_o);
out.H_y = tseries(dateFrom:dateTo,H_y);
out.D_o = tseries(dateFrom:dateTo,D_o);
out.D_y = tseries(dateFrom:dateTo,D_y);
out.X = out.X_o+out.X_y;
out.I = out.I_o+out.I_y;
out.H = out.H_o+out.H_y;
out.D = out.D_o+out.D_y;

end