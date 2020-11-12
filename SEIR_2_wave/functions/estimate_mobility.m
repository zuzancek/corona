function []=estimate_mobility(mob,Rt,dateFrom,delay)

figure;
plot(mob/100,'linewidth',1); hold on;
plot(Rt,'linewidth',1);
grid on;
legend({'mobility','Rt'});

end