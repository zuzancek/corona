fuunction []=simulate_model(

s = setparam();
vacc_plan;


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
    E_o(t+1,:)=E_o(t,:)+S_o(t,:).*mu_o.*R0./T_inf.*F(t)-E_o(t,:)./T_lat;
    E_y(t+1,:)=E_y(t,:)+S_y(t,:).*mu_y.*R0./T_inf.*F(t)-E_y(t,:)./T_lat;
end