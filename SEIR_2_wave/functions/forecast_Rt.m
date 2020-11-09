function [Rt,q_mat,It,Xt,x_mat] = forecast_Rt(Rt_old, dI_old,per, I0,pop_size,T_rem,N,q_vec)

figure;
plot(dI_old);hold on;


figure;
plot(Rt_old);hold on;

n = length(dI_old);
ddI = dI_old(2:end)-dI_old(1:end-1);
figure;
plot(ddI);hold on;
% % x = 1:n-1;
% % xx = 1:n-1+per;
% % y = interp1(x,ddI,xx,'spline');
% % plot(y);


end