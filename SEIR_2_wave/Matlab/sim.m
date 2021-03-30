function [] = simulate_model(init,vacc_plan)

s = setparam();
% vacc_plan;

k_inf = s.k_inf;
gamma_lat = repmat((s.pdf_lat(end:-1:1)/(k_inf:-1:0))',1,3);
k_pre = s.k_pre; T_test = s.T_test.mean;
sigma_o = s.obs_ratio; sigma_y = s.obs_ratio;
gamma_test = repmat((s.pdf_pre(end:-1:1)/((k_pre:-1:0)+T_test))',1,3);
k_inf = s.k_inf;
gamma_inf = repmat((s.pdf_inf(end:-1:1)/(k_inf:-1:0))',1,3);

for t=t0+1:T
    % vaccination
    vacc_sus_1_o = S_o(1,t)/(S_o(1,t)+R_o(1,t));
    vacc_sus_2_o = S_o(2,t)/(S_o(2,t)+R_o(2,t));
    vaccine_intake_sus_o(t,:) = [-vacc_sus_1_o.*vacc_plan.D1_o;vacc_sus_1_o.*vacc_plan.D1_o-vacc_sus_2_o.*vacc_plan.D2_o;vacc_sus_2_o.*vacc_plan.D2_o];
    vaccine_intake_rem_o(t,:) = [-(1-vacc_sus_1_o).*vacc_plan.D1_o;(1-vacc_sus_1_o).*vacc_plan.D1_o-(1-vacc_sus_2_o).*vacc_plan.D2_o;(1-vacc_sus_2_o).*vacc_plan.D2_o];
    vacc_sus_1_y = S_y(1,t)/(S_y(1,t)+R_y(1,t));
    vacc_sus_2_y = S_y(2,t)/(S_y(2,t)+R_y(2,t));
    vaccine_intake_sus_y(t,:) = [-vacc_sus_1_y.*vacc_plan.D1_y;vacc_sus_1_y.*vacc_plan.D1_y-vacc_sus_2_y.*vacc_plan.D2_y;vacc_sus_2_y.*vacc_plan.D2_y];
    vaccine_intake_rem_y(t,:) = [-(1-vacc_sus_1_y).*vacc_plan.D1_y;(1-vacc_sus_1_y).*vacc_plan.D1_y-(1-vacc_sus_2_y).*vacc_plan.D2_y;(1-vacc_sus_2_y).*vacc_plan.D2_y];
    % suspectible
    S_o(t+1,:)=S_o(t,:)+vaccine_intake_sus_o(t,:)-S_o(t,:).*mu_o.*R0./T_inf.*F(t);
    S_y(t+1,:)=S_y(t,:)+vaccine_intake_sus_y(t,:)-S_y(t,:).*mu_y.*R0./T_inf.*F(t);
    % exposed    
    E_o(t+1,:)=E_o(t,:)+S_o(t,:).*mu_o.*R0./T_inf.*F(t)-dot(E_o(t-k_inf:t-1,:),gamma_lat);
    E_y(t+1,:)=E_y(t,:)+S_y(t,:).*mu_y.*R0./T_inf.*F(t)-dot(E_y(t-k_inf:t-1,:),gamma_lat);
    % infectious, asymptomatic/mild - unobserved
    X_o(t,:) = sigma_o*dot(A_o(t-k_pre-T_test:t-1-T_test,:),gamma_test);
    X_y(t,:) = sigma_y*dot(A_y(t-k_pre-T_test:t-1-T_test,:),gamma_test);
    A_o(t+1,:)=A_o(t,:)+dot(E_o(t-k_inf:t-1,:),gamma_lat)-X_o(t,:)-(1-sigma_o)*dot(A_o(t-k_inf:t-1,:),gamma_inf);
    A_y(t+1,:)=A_y(t,:)+dot(E_y(t-k_inf:t-1,:),gamma_lat)-X_y(t,:)-(1-sigma_y)*dot(A_y(t-k_inf:t-1,:),gamma_inf);
    % infectious (at least at the beginning), observed
    I_o(t+1,:)=I_o(t,:)+X_o(t,:)-eta_o*dot(I_o(t-k_hosp:t-1,:),gamma_test)-(1-eta_o)*dot(A_o(t-k_inf:t-1,:),gamma_inf);
    I_y(t+1,:)=I_y(t,:)+X_y(t,:)-eta_y*dot(I_y(t-k_hosp:t-1,:),gamma_test)-(1-eta_y)*dot(A_y(t-k_inf:t-1,:),gamma_inf);
    
end

end