function [x] = get_first_infections(obs,obs_ratio,N)

total_number = min(N,obs/obs_ratio);
x = zeros(N,1);
x(1:total_number) = 1;
x = x(randperm(N));

end

