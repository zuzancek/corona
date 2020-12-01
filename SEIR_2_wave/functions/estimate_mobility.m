function [res] = estimate_mobility(mob,Rt,dateFrom,delay,dateFrom_fullSample)

r_dateFrom = dateFrom;
r_dateTo = enddate(Rt);
r = resize(Rt,r_dateFrom:r_dateTo);
r_full = double(resize(Rt,dateFrom_fullSample:r_dateTo));
r_data = double(r);

m_dateFrom = dateFrom-delay;
m_dateTo = r_dateTo-delay;
m = resize(mob.orig/100,m_dateFrom:m_dateTo);
% m_raw = double(resize(mob.raw/100,m_dateFrom:m_dateTo));
m_full = double(resize(mob.orig/100,dateFrom_fullSample-delay:m_dateTo));
m_data = double(m);

[~,m_dmin] = min(m);
[~,r_dmin] = min(r);
[r_vmax,r_dmax] = max(r_full);
% assert(r_dmin-m_dmin==delay);

figure('Name','Mobility as a predictor of Rt');
plot(m,'linewidth',2); hold on;
plot(r,'linewidth',2);
grid on;
legend({'mobility','Rt'});

% curve fitting: R(m) = (1/(a_0-m^k)-a_1)*c_2
m_idx = find(m_data<1,1); a_m = interp1([m_data(m_idx-1) m_data(m_idx)],[0 1],1);
r_1 = r_data(m_idx-1)*a_m+(1-a_m)*r_data(m_idx);
r_idx = find(r_data<1,1); 
if isempty(r_idx)
    [~,r_idx] = min(r_data);
end
a_r = interp1([r_data(r_idx-1) r_data(r_idx)],[0 1],1);
m_1 = m_data(r_idx-1)*a_r+(1-a_r)*m_data(r_idx);
m_2 = m_full(r_dmax);
sl = (r_vmax-r_1)/(m_2-1);

figure('Name','Mobility, impact on Rt (raw data)');
x1 = 0.65:0.01:1;
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
grid on;
xlabel('Mobility (%)');
ylabel('R_t');

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
        a = exp((m_1-1)*k)*(r_1-1)/(1-r_1*exp((m_1-1)*k));
        yy = r_1*(a+1).*(1-a./(a+exp((x-1).*k)));
    end

end