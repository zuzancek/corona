function [res_def,res_real] = estimate_mobility(mob,inputs,dateFrom,delay,dateFrom_fullSample)

delay = 15;
delay_real = 10;
Rt = inputs.default;
Rt_real = inputs.real;
r_dateFrom = dateFrom;
r_dateTo = enddate(Rt);
r = resize(Rt,r_dateFrom:r_dateTo);
r_real = resize(Rt_real,r_dateFrom:r_dateTo);
r_full = double(resize(Rt,dateFrom_fullSample:r_dateTo));
r_data = double(r);
r_real_full = double(resize(Rt_real,dateFrom_fullSample:r_dateTo));
r_real_data = double(r_real);

m_dateFrom = dateFrom-delay;
m_dateTo = r_dateTo-delay;
m = resize(mob.orig/100,m_dateFrom:m_dateTo);
% m_raw = double(resize(mob.raw/100,m_dateFrom:m_dateTo));
m_full = double(resize(mob.orig/100,dateFrom_fullSample-delay:m_dateTo));
m_data = double(m);

m_real_dateFrom = dateFrom-delay_real;
m_real_dateTo = r_dateTo-delay_real;
m_real = resize(mob.orig/100,m_real_dateFrom:m_real_dateTo);
% m_raw = double(resize(mob.raw/100,m_dateFrom:m_dateTo));
m_real_full = double(resize(mob.orig/100,dateFrom_fullSample-delay_real:m_real_dateTo));
m_real_data = double(m_real);

[m_v_min,m_dmin] = min(m);
[r_v_min,r_dmin] = min(r);
[m_real_v_min,m_real_dmin] = min(m_real);
[r_real_v_min,r_real_dmin] = min(r_real);
[r_real_vmax,r_real_dmax] = max(r_real_full);
fprintf('Delay in reported data: %d\n',r_dmin-m_dmin);
fprintf('Delay in real data: %d\n',r_real_dmin-m_real_dmin);
% assert(r_dmin-m_dmin==delay);

figure('Name','Mobility as a predictor of Rt');
plot(m,'linewidth',2); hold on;
plot(r,'linewidth',2);
plot(r_real,'linewidth',2);
grid on;
legend({'mobility','Rt','Rt_{real}'});

% curve fitting: R(m) = (1/(a_0-m^k)-a_1)*c_2
m_idx = find(m_data<1,1); a_m = interp1([m_data(m_idx-1) m_data(m_idx)],[0 1],1);
r_1 = r_data(m_idx-1)*a_m+(1-a_m)*r_data(m_idx);
r_real_1 = r_real_data(m_idx-1)*a_m+(1-a_m)*r_real_data(m_idx);
m_real_2 = m_full(r_real_dmax);
sl = (r_real_vmax-r_real_1)/(m_real_2-1);

figure('Name','Mobility, impact on Rt');
pp1=scatter(m_full,r_full); hold on;
pp2=scatter(m_full,r_real_full); hold on;
x0_axis = m_v_min:0.01:1.01;
nn = length(x0_axis);
k0_pos = -0.45; k0_real_pos = -0.65;
k0_neg = 3; k0_real_neg = 6;
[y_def_pos]=R_fun(m_v_min,r_v_min,r_1,x0_axis(1:end),k0_pos);
[y_def_neg]=R_fun_neg(m_v_min,r_v_min,r_1,x0_axis(1:end),k0_neg);
plot(x0_axis(1:end),y_def_pos,'-','linewidth',2,'Color',pp1.CData,'Marker','<','MarkerIndices',[1:12:nn]); %#ok<*NBRAK>
plot(x0_axis(1:end),y_def_neg,'--','linewidth',2,'Color',pp1.CData,'Marker','>','MarkerIndices',[1:12:nn]);
[y_real_pos]=R_fun(m_v_min,r_real_v_min,r_real_1,x0_axis(1:end),k0_real_pos);
[y_real_neg]=R_fun_neg(m_v_min,r_real_v_min,r_real_1,x0_axis(1:end),k0_real_neg);
plot(x0_axis(1:end),y_real_pos,'-','linewidth',2,'Color',pp2.CData,'Marker','<','MarkerIndices',[1:12:nn]);
plot(x0_axis(1:end),y_real_neg,'--','linewidth',2,'Color',pp2.CData,'Marker','>','MarkerIndices',[1:12:nn]);
x1_axis = 1:0.01:1.45;
y2 = r_1+sl.*(x1_axis-x1_axis(1));
plot(x1_axis,y2+0.001,'Color',pp1.CData,'linewidth',2);
plot(x1_axis,y2,'Color',pp2.CData,'linewidth',2);
grid on;
xlabel('mobility');
ylabel('Rt')
legend([pp1,pp2],{'reported','real'});
title('Raw data');

figure('Name','Mobility, impact on Rt (smooth)');
x = [x0_axis,x1_axis];
subplot(2,1,1)
pp1=plot(x,[y_def_pos,y2+0.01],'linewidth',0.5,'linestyle',':');hold on;
pp2=plot(x,[y_def_neg,y2],'linewidth',0.5,'linestyle',':');hold on;
p1=plot(x,smooth_series([y_def_pos,y2],5,5,1)+0.01,'linewidth',1,'Color',pp1.Color,'Marker','<','MarkerIndices',[1:10:nn]);hold on;
p2=plot(x,smooth_series([y_def_neg,y2],5,5,1),'linewidth',1,'Color',pp2.Color,'Marker','>','MarkerIndices',[1:10:nn]);hold on;
grid on;
xlabel('mobility');
ylabel('Rt')
legend([p1,p2],{'increasing mobility','decreasing mobility'});
title('Reported data');

subplot(2,1,2)
pp1=plot(x,[y_real_pos,y2+0.01],'linewidth',0.5,'linestyle',':');hold on;
pp2=plot(x,[y_real_neg,y2],'linewidth',0.5,'linestyle',':');hold on;
p1=plot(x,smooth_series([y_real_pos,y2+0.01],5,5,1),'linewidth',1,'Color',pp1.Color,'Marker','<','MarkerIndices',[1:10:nn]);hold on;
p2=plot(x,smooth_series([y_real_neg,y2],5,5,1),'linewidth',1,'Color',pp2.Color,'Marker','>','MarkerIndices',[1:10:nn]);hold on;
grid on;
xlabel('mobility');
ylabel('Rt')
legend([p1,p2],{'increasing mobility','decreasing mobility'});
title('Real (model-implied) data');

res_def = struct;
res_def.x = x;
res_def.y_pos = smooth_series([y_def_pos,y2],5,5,1);
res_def.y_neg = smooth_series([y_def_neg,y2],5,5,1);
res_real = struct;
res_real.x = x;
res_real.y_pos = smooth_series([y_real_pos,y2],5,5,1);
res_real.y_neg = smooth_series([y_real_neg,y2],5,5,1);

    function [yy,aa0]=R_fun(xx0,yy0,yy1,x,k0)
        aa0 = (yy1-yy0)/(1-xx0.^k0);
        bb0 = yy1-aa0;
        yy = aa0.*x.^k0+bb0;
    end
    function [yy,aa0]=R_fun_neg(xx0,yy0,yy1,x,k0)
        aa0 = (yy1-yy0)/(1-exp((xx0-1)*k0));
        bb0 = yy1-aa0;
        yy = aa0.*exp(k0*(x-1))+bb0;
    end
end