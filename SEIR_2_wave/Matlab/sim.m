function [] = simulate_model(init,vacc_plan)

s = setparam();
% vacc_plan;

k_inf = s.k_inf;
gamma_ea = repmat((s.pdf_lat(end:-1:1)/(k_inf:-1:0))',1,3);
k_pre = s.k_pre; T_test = s.T_test.mean;
sigma_o = s.obs_ratio; sigma_y = s.obs_ratio;
gamma_ai = repmat((s.pdf_pre(end:-1:1)/((k_pre:-1:0)+T_test))',1,3);
k_inf = s.k_inf;
gamma_ar = repmat((s.pdf_inf(end:-1:1)/(k_inf:-1:0))',1,3);
k_hosp = s.k_hosp; 
eta_o = s.eta_o; eta_y = s.eta_y;
gamma_iho = repmat((s.pdf_ih_o(end:-1:1)/(k_hosp:-1:0))',1,3);
gamma_ihy = repmat((s.pdf_ih_y(end:-1:1)/(k_hosp:-1:0))',1,3);
k_reci = s.k_sick;
gamma_iro = repmat((s.pdf_ir_o(end:-1:1)/(k_reci:-1:0))',1,3);
gamma_iry = repmat((s.pdf_ir_y(end:-1:1)/(k_reci:-1:0))',1,3);
k_death = s.k_death;
omega_o = s.omega_o;omega_y = s.omega_y;
gamma_hdo = repmat((s.pdf_hd_o(end:-1:1)/(k_death:-1:0))',1,3);
gamma_hdy = repmat((s.pdf_hd_y(end:-1:1)/(k_death:-1:0))',1,3);
k_rec = s.k_rec; 
gamma_hro = repmat((s.pdf_hr_o(end:-1:1)/(k_rec:-1:0))',1,3);
gamma_hry = repmat((s.pdf_hr_y(end:-1:1)/(k_rec:-1:0))',1,3);

% **** arrays
S_o = zeros(T_total,3); S_y = S_o; E_y = S_o; E_o = S_o; 
X_o = S_o; X_y = S_o; A_o = S_o; A_y = S_o; I_o = S_o; I_y = S_o;
D_o = S_o; D_y = S_o; H_o = S_o; H_y = S_o; T_o = S_o; T_y = S_o; 
R_o = S_o; R_y = S_o; P_o = S_o; P_y = S_o;
d_I_H_o = S_o; d_I_H_y = S_o; d_H_D_o = S_o; d_H_D_y = S_o;
d_H_R_o = S_o; d_H_R_y = S_o;
vaccine_intake_sus_o = S_o; vaccine_intake_rem_o = S_o;
vaccine_intake_sus_y = S_o; vaccine_intake_rem_y = S_o;

for t=t0+1:T_total
    % vaccination
    vacc_sus_1_o = S_o(1,t)/(S_o(1,t)+T_o(1,t));
    vacc_sus_2_o = S_o(2,t)/(S_o(2,t)+T_o(2,t));
    vaccine_intake_sus_o(t-1,:) = [-vacc_sus_1_o.*vacc_plan.D1_o;vacc_sus_1_o.*vacc_plan.D1_o-vacc_sus_2_o.*vacc_plan.D2_o;vacc_sus_2_o.*vacc_plan.D2_o];
    vaccine_intake_rem_o(t-1,:) = [-(1-vacc_sus_1_o).*vacc_plan.D1_o;(1-vacc_sus_1_o).*vacc_plan.D1_o-(1-vacc_sus_2_o).*vacc_plan.D2_o;(1-vacc_sus_2_o).*vacc_plan.D2_o];
    vacc_sus_1_y = S_y(1,t)/(S_y(1,t)+T_y(1,t));
    vacc_sus_2_y = S_y(2,t)/(S_y(2,t)+T_y(2,t));
    vaccine_intake_sus_y(t-1,:) = [-vacc_sus_1_y.*vacc_plan.D1_y;vacc_sus_1_y.*vacc_plan.D1_y-vacc_sus_2_y.*vacc_plan.D2_y;vacc_sus_2_y.*vacc_plan.D2_y];
    vaccine_intake_rem_y(t-1,:) = [-(1-vacc_sus_1_y).*vacc_plan.D1_y;(1-vacc_sus_1_y).*vacc_plan.D1_y-(1-vacc_sus_2_y).*vacc_plan.D2_y;(1-vacc_sus_2_y).*vacc_plan.D2_y];
    % suspectible
    d_S_E_o=S_o(t-1,:).*mu_o.*R0./T_inf.*F(t);
    d_S_E_y=S_y(t-1,:).*mu_y.*R0./T_inf.*F(t);
    S_o(t,:) = S_o(t-1,:)+vaccine_intake_sus_o(t-1,:)-d_S_E_o;
    S_y(t,:) = S_y(t-1,:)+vaccine_intake_sus_y(t-1,:)-d_S_E_y;
    % exposed    
    d_E_A_o=dot(E_o(t-k_inf:t-1,:),gamma_ea);
    d_E_A_y=dot(E_y(t-k_inf:t-1,:),gamma_ea);
    E_o(t,:) = E_o(t-1,:)+d_S_E_o-d_E_A_o;
    E_y(t,:) = E_y(t-1,:)+d_S_E_y-d_E_A_y;
    % infectious, asymptomatic/mild - unobserved
    X_o(t,:) = sigma_o*dot(A_o(t-k_pre-T_test:t-1-T_test-1,:),gamma_ai);
    X_y(t,:) = sigma_y*dot(A_y(t-k_pre-T_test:t-1-T_test-1,:),gamma_ai);
    d_A_R_o=(1-sigma_o)*dot(A_o(t-k_inf:t-1,:),gamma_ar);
    d_A_R_y=(1-sigma_y)*dot(A_y(t-k_inf:t-1,:),gamma_ar);
    A_o(t,:) = A_o(t-1,:)+d_E_A_o-X_o(t,:)-d_A_R_o;
    A_y(t,:) = A_y(t-1,:)+d_E_A_y-X_y(t,:)-d_A_R_y;
    % infectious (at least at the beginning), observed
    d_I_H_o(t,:) = eta_o*dot(I_o(t-k_hosp:t-1,:),gamma_iho);
    d_I_H_y(t,:) = eta_y*dot(I_y(t-k_hosp:t-1,:),gamma_ihy);
    d_I_R_o = (1-eta_o)*dot(I_o(t-k_reci:t-1,:),gamma_iro);
    d_I_R_y = (1-eta_y)*dot(I_y(t-k_reci:t-1,:),gamma_iry);
    I_o(t,:) = I_o(t-1,:)+X_o(t-1,:)-d_I_H_o(t,:)-d_I_R_o;
    I_y(t,:) = I_y(t-1,:)+X_y(t-1,:)-d_I_H_y(t,:)-d_I_R_y;    
    % hospitalisations
    d_H_D_o(t,:) = omega_o*dot(H_o(t-k_death:t-1,:),gamma_hdo);
    d_H_D_y(t,:) = omega_y*dot(H_y(t-k_death:t-1,:),gamma_hdy);
    d_H_R_o(t,:) = (1-omega_o)*dot(H_o(t-k_death:t-1,:),gamma_hdo);
    d_H_R_y(t,:) = (1-omega_y)*dot(H_y(t-k_death:t-1,:),gamma_hdy);
    H_o(t,:) = H_o(t-1,:)+d_I_H_o(t,:)-d_H_D_o(t,:)-d_H_R_o(t,:);
    H_y(t,:) = H_y(t-1,:)+d_I_H_y(t,:)-d_H_D_y(t,:)-d_H_R_y(t,:);
    % deaths
    D_o(t,:) = D_o(t-1,:)+d_H_D_o(t,:);
    D_y(t,:) = D_y(t-1,:)+d_H_D_y(t,:);
    % (temporarily/permanently) removed
    T_o(t,:) = T_o(t-1,:)+vaccine_intake_rem_o(t-1,:).*T_o(t-1,:)./R_o(t-1,:)+d_I_R_o+d_A_R_o;
    T_y(t,:) = T_y(t-1,:)+vaccine_intake_rem_y(t-1,:).*T_y(t-1,:)./R_y(t-1,:)+d_I_R_y+d_A_R_y;
    P_o(t,:) = P_o(t-1,:)+vaccine_intake_rem_o(t-1,:).*P_o(t-1,:)./R_o(t-1,:)+d_H_R_o(t,:);
    P_y(t,:) = P_y(t-1,:)+vaccine_intake_rem_y(t-1,:).*P_y(t-1,:)./R_y(t-1,:)+d_H_R_y(t,:);
    R_o(t,:) = P_o(t,:)+T_o(t,:);R_y(t,:) = P_y(t,:)+T_y(t,:);
end

end