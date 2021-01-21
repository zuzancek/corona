function [out] = XIHDt(x,p,s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
delay = max([s.k_hosp,s.k_sick,s.k_death,s.k_rec])+1;
T_total = T+delay;

I0 = init.I;     
H0 = init.H; 
D0 = init.D;

method_params = s.smoothing_method_params;
try
    rho = init.rho(end-T_total+1:end);
catch err
    rho = method_params(init.rho); 
    rho = rho(dateFrom-delay:dateTo);
end
try
    T_delay = p.T_delay(end-T_total+1:end);
catch err
    T_delay = zeros(T_total,1);
end
varsigma = method_params(init.varsigma); varsigma = varsigma(dateFrom:dateTo);
     
% ********* initialization
method_data = s.smoothing_method_data;
X = method_data(x.NewCases(dateFrom:dateTo));

% ********* arrays
X_o = X.*rho(delay+1:end);   X_y = X-X_o;
D0 = method_data(D0(dateFrom-delay+1:dateTo));
D_o = zeros(T_total,1); D_o(1:delay) = D0(1:delay)*varsigma(1); 
D_y = zeros(T_total,1); D_y(1:delay) = D0(1:delay)-D_o(1:delay); 
H0 = method_data(H0(dateFrom-delay:dateTo));
H_o = zeros(T_total,1); H_o(1:delay) = H0(1:delay).*s.eta_o.*rho(1:delay)./(s.eta_o.*rho(1:delay)+s.eta_y.*(1-rho(1:delay)));
H_y = zeros(T_total,1); H_y(1:delay) = H0(1:delay)-H_o(1:delay);   
I0 = method_data(I0(dateFrom-delay+1:dateTo));
I_o = zeros(T_total,1); I_o(1:delay) = I0(1:delay).*rho(1:delay);
I_y = zeros(T_total,1); I_y(1:delay) = I0(1:delay)-I_o(1:delay);

d_I_H_o = zeros(T_total,1);  d_I_H_y = d_I_H_o; d_I_R_o = d_I_H_o; d_I_R_y = d_I_H_o;
d_H_D_o = d_I_H_o; d_H_D_y = d_I_H_o; d_H_R_o = d_I_H_o; d_H_R_y = d_I_H_o;

% ******* parameters
% hospital admission
k_hosp = s.k_hosp;      x_hosp = (1:k_hosp)'; 
T_hosp_o = s.T_hosp_o;  p_T_hosp_o = pdf(s.T_hosp_pdf_type,x_hosp,T_hosp_o+zeros(k_hosp,1));
T_hosp_y = s.T_hosp_y;  p_T_hosp_y = pdf(s.T_hosp_pdf_type,x_hosp,T_hosp_y+zeros(k_hosp,1));
eta_o = s.eta_o;        alpha_iho = eta_o./x_hosp.*p_T_hosp_o;  alpha_iho = alpha_iho(end:-1:1);
eta_y = s.eta_y;        alpha_ihy = eta_y./x_hosp.*p_T_hosp_y;  alpha_ihy = alpha_ihy(end:-1:1);
% recovery from sickness at home
k_sick = s.k_sick;      x_sick = (1:k_sick)';       x_sick_mat = repmat(x_sick',T_total,1); T_sick_std2 = s.T_sick_std^2;
T_sick_o = s.T_sick_o-T_delay-p.T_test_to_result;   T_sick_o_shape = T_sick_o*T_sick_std2; T_sick_scale = 1/T_sick_std2;
p_T_sick_o = pdf(s.T_sick_pdf_type,x_sick_mat,repmat(T_sick_o_shape,1,k_sick),repmat(T_sick_scale,T_total,k_sick));
T_sick_y = s.T_sick_y-T_delay-p.T_test_to_result;   T_sick_y_shape = T_sick_y*T_sick_std2; 
p_T_sick_y = pdf(s.T_sick_pdf_type,x_sick_mat,repmat(T_sick_y_shape,1,k_sick),repmat(T_sick_scale,T_total,k_sick));
alpha_iro = (1-eta_o)./x_sick_mat.*p_T_sick_o;      alpha_iro = alpha_iro(:,end:-1:1);
alpha_iry = (1-eta_y)./x_sick_mat.*p_T_sick_y;      alpha_iry = alpha_iry(:,end:-1:1);
% death at hospital
k_death = s.k_death;        x_death = (1:k_death)';
T_death_o = s.T_death_o;    p_T_death_o = pdf(s.T_death_pdf_type,x_death,T_death_o+zeros(k_death,1));
T_death_y = s.T_death_y;    p_T_death_y = pdf(s.T_death_pdf_type,x_death,T_death_y+zeros(k_death,1));
omega_o = s.omega_o;        alpha_hdo = omega_o./x_death.*p_T_death_o;  alpha_hdo = alpha_hdo(end:-1:1);
omega_y = s.omega_y;        alpha_hdy = omega_y./x_death.*p_T_death_y;  alpha_hdy = alpha_hdy(end:-1:1);
% recovery at hospital
k_rec = s.k_rec;            x_rec = (1:k_rec)';                     T_rec_std2 = s.T_rec_std^2; T_rec_scale = 1/T_rec_std2+zeros(k_rec,1);
T_rec_o = s.T_rec_o;        T_rec_o_shape = T_rec_o*T_sick_std2;    
p_T_rec_o = pdf(s.T_rec_pdf_type,x_rec,T_rec_o_shape+zeros(k_rec,1),T_rec_scale);
T_rec_y = s.T_rec_y;        T_rec_y_shape = T_rec_y*T_sick_std2;    
p_T_rec_y = pdf(s.T_rec_pdf_type,x_rec,T_rec_y_shape+zeros(k_rec,1),T_rec_scale);
alpha_hro = (1-omega_o)./x_rec.*p_T_rec_o;          alpha_hro = alpha_hro(end:-1:1);
alpha_hry = (1-omega_y)./x_rec.*p_T_rec_y;          alpha_hry = alpha_hry(end:-1:1);

% ******* Equations
% I(t+1) = I(t)+X(t)-I_H(t)-I_R(t);     
% H(t+1) = H(t)+I_H(t)-H_D(t)-H_R(t);
% D(t+1) = D(t)+H_D(t);

% ******* Calculation
for t=delay:T_total-1 
    d_I_H_o(t) = dot(alpha_iho,I_o(t-k_hosp+1:t));     d_I_R_o(t) = dot(alpha_iro(t,:),I_o(t-k_sick+1:t));
    I_o(t+1) = I_o(t)+X_o(t-delay+1)-d_I_H_o(t)-d_I_R_o(t);
    d_I_H_y(t) = dot(alpha_ihy,I_y(t-k_hosp+1:t));     d_I_R_y(t) = dot(alpha_iry(t,:),I_y(t-k_sick+1:t));
    I_y(t+1) = I_y(t)+X_y(t-delay+1)-d_I_H_y(t)-d_I_R_y(t);
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
out.I_o = tseries(dateFrom:dateTo,I_o(delay+1:end));
out.I_y = tseries(dateFrom:dateTo,I_y(delay+1:end));
out.H_o = tseries(dateFrom:dateTo,H_o(delay+1:end));
out.H_y = tseries(dateFrom:dateTo,H_y(delay+1:end));
out.D_o = tseries(dateFrom:dateTo,D_o(delay+1:end));
out.D_y = tseries(dateFrom:dateTo,D_y(delay+1:end));
out.X = out.X_o+out.X_y;
out.I = out.I_o+out.I_y;
out.H = out.H_o+out.H_y;
out.D = out.D_o+out.D_y;

end