function [res] = estimate_mobility(mob,inputs,dI,dateFrom,delay,dateFrom_fullSample)

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

di = resize(dI,r_dateFrom:r_dateTo); 
di = di/di(r_dateFrom);

[m_v_min,m_dmin] = min(m);
[r_v_min,r_dmin] = min(r);
[r_real_v_min,r_real_dmin] = min(r_real);
[r_vmax,r_dmax] = max(r_full);
[r_real_vmax,r_real_dmax] = max(r_real_full);
fprintf('Delay in reported data: %d\n',r_dmin-m_dmin);
fprintf('Delay in real data: %d\n',r_real_dmin-m_dmin);
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
r_idx = find(r_data<1,1); 
if isempty(r_idx)
    [~,r_idx] = min(r_data);
end
r_real_1 = r_real_data(m_idx-1)*a_m+(1-a_m)*r_real_data(m_idx);
r_real_idx = find(r_real_data<1,1); 
if isempty(r_real_idx)
    [~,r_real_idx] = min(r_real_data);
end
a_r = interp1([r_data(r_idx-1) r_data(r_idx)],[0 1],1);
m_1 = m_data(r_idx-1)*a_r+(1-a_r)*m_data(r_idx);
m_2 = m_full(r_dmax);
sl = (r_vmax-r_1)/(m_2-1);
a_r_real = interp1([r_real_data(r_real_idx-1) r_real_data(r_real_idx)],[0 1],1);
m_real_1 = m_data(r_real_idx-1)*a_r_real+(1-a_r_real)*m_data(r_real_idx);
m_real_2 = m_full(r_real_dmax);
sl = (r_real_vmax-r_real_1)/(m_real_2-1);

figure('Name','Mobility, impact on Rt (raw data)');
pp1=scatter(m_full,r_full); hold on;
pp2=scatter(m_full,r_real_full); hold on;

x1 = 0.65:0.01:1;
x_axis = m_v_min:0.01:1;
k0_pos = -0.45; k0_real_pos = -0.65;
k0_neg = 3; k0_real_neg = 6;
[y_def_pos]=R_fun(m_v_min,r_v_min,r_1,x_axis(2:end),k0_pos);
[y_def_neg]=R_fun_neg(m_v_min,r_v_min,r_1,x_axis(2:end),k0_neg);
plot(x_axis(2:end),y_def_pos,'-','linewidth',2,'Color',pp1.CData);
plot(x_axis(2:end),y_def_neg,'--','linewidth',2,'Color',pp1.CData);
[y_real_pos]=R_fun(m_v_min,r_real_v_min,r_real_1,x_axis(2:end),k0_real_pos);
[y_real_neg]=R_fun_neg(m_v_min,r_real_v_min,r_real_1,x_axis(2:end),k0_real_neg);
plot(x_axis(2:end),y_real_pos,'-','linewidth',2,'Color',pp2.CData);
plot(x_axis(2:end),y_real_neg,'--','linewidth',2,'Color',pp2.CData);
grid on;
legend([pp1,pp2],{'reported','real'});


x2 = 1:0.01:1.45;
y2 = r_1+sl.*(x2-x2(1));

k_array = [0.012:0.2:2];
idx_chosen = 1;
for i=1:length(k_array)
    k = k_array(i);
    y1 = modelfun(k,x1);
    if i~=idx_chosen
        plot(x1,y1,'linewidth',0.5,'Color',[1-i*0.075 0.5 i*0.075]);hold on;
    end
end
k = k_array(idx_chosen);
y1 = modelfun2(k,x1);
plot([x1,x2],[y1,y2],'linewidth',2,'Color','m');hold on;
scatter(m_full,r_full,'k'); hold on;
scatter(m_full,r_real_full,'g'); hold on;
grid on;
xlabel('Mobility (%)');
ylabel('R_t');

%
[~,p0] = min(m_full);
xm_1 = m_full(1:p0); yr_1 = r_full(1:p0);
xm_2 = m_full(p0:end); yr_2 = r_full(p0:end);
iidx = find(m_full(1:p0)<1,1);
xm_2(end+1) = 1; yr_2(end+1) = r_full(iidx);
yr_3 = interp1(xm_1,yr_1,xm_2);
yr_mean = 0.5*(yr_2+yr_3);
figure;
plot(xm_1,yr_1);hold on;
plot(xm_2,yr_2);hold on;
plot(xm_2,yr_mean);hold on;
grid on;

figure('Name','Mobility, impact on Rt (smooth function)');
x = [x1 x2(2:end)];
y = [y1 y2(2:end)];
plot(x,y,'linewidth',1); hold on;
z = smooth_series(y,5,5,1);
plot(x,z,'linewidth',2); hold on;
grid on;
legend({'raw','smooth'})
xlabel('Mobility (%)');
ylabel('R_t');

res = struct();
res.f1.x = x1;
[yy,a] = modelfun(k_array(idx_chosen),x1);
res.f1.y = yy;
res.f1.p = [r_1 a k_array(idx_chosen)];
res.f2.x = x2;
res.f2.y = y2;
res.f2.s = sl;
res.f2.k = r_1;
res.f.x_grid = x;
res.f.y_grid = z;

    function [yy,a]=modelfun(k,x)
        a = m_1^k*(r_1-1)/(1-r_1*m_1^k);
        yy = r_1*(a+1).*(1-a./(a+x.^k));
    end
    function [yy,a]=modelfun2(k,x)
        a = exp((m_1^3-1)*k)*(r_1-1)/(1-r_1*exp((m_1^3-1)*k));
        yy = r_1*(a+1).*(1-a./(a+exp((x.^3-1).*k)));
    end

%     function [yy,aa0]=R_fun0(xx0,yy0,yy1,x,k0)
%         aa0 = fzero(@inner_fnc,0.5);
%         function [d]=inner_fnc(aa)
%             d=aa./(aa-1).*yy0/yy1-(1./(aa-xx0.^k0)-1);
%         end
%         yy = yy1*(aa0-1)./aa0.*(1./(aa0-x.^k0)-1);
%     end
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