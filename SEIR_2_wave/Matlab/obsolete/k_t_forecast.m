

function [fcast] = k_t_forecast(eta,omega,phi) 
% k, k_t_star_pred, k_t_constant_dow
lambda = 1+k*(phi-1)/30;
w = zeros(length(k),1);
w(k<=omega+1) = 1-((k(k<=omega+1)-1)/omega).^2;
eta_star = median_k_t.*eta;
x = [eta_star,k_t_star_pred];
fcast = lambda.*(w.*min(1,x)+(1-w).*k_t_constant_dow);
end

function [dist] = joint_distr(eta,omega,phi)
    k_t_forecast_val = k_t_forecast(eta,omega,phi);
    
