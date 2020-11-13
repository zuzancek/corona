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
assert(r_dmin-m_dmin==delay);

%
figure;
plot(m,'linewidth',1); hold on;
plot(r,'linewidth',1);
grid on;
legend({'mobility','Rt'});

% curve fitting: R(m) = (1/(a_0-m^k)-a_1)*c_2
m_idx = find(m_data<1,1); a_m = interp1([m_data(m_idx-1) m_data(m_idx)],[0 1],1);
r_1 = r_data(m_idx-1)*a_m+(1-a_m)*r_data(m_idx);
r_idx = find(r_data<1,1); a_r = interp1([r_data(r_idx-1) r_data(r_idx)],[0 1],1);
m_1 = m_data(r_idx-1)*a_r+(1-a_r)*m_data(r_idx);
m_2 = m_full(r_dmax);
sl = (r_vmax-r_1)/(m_2-1);

figure;
x1 = 0.65:0.01:1;
x2 = 1:0.01:1.45;
y2 = r_1+sl.*(x2-x2(1));

k_array = [1.5:0.5:2.5];
for i=1:length(k_array)
    k = k_array(i);
    y1 = modelfun(k,x1);
    plot(x1,y1,'linewidth',1);hold on;
end
plot(x2,y2,'linewidth',1);hold on;
scatter(m_full,r_full,'k'); hold on;
grid on;
xlabel('Mobility (%)');
ylabel('R_t');

figure;
x = [x1 x2(2:end)];
y = [y1 y2(2:end)];
plot(x,y,'linewidth',1); hold on;
z = smooth_series(y,7,5,1);
plot(x,z,'linewidth',1); hold on;
grid on;
xlabel('Mobility (%)');
ylabel('R_t');

res = struct();
res.f1.x = x1;
[yy,a] = modelfun(k_array(1),x1);
res.f1.y = yy;
res.f1.p = [r1 a k_array(1)];
res.f2.x = x2;
res.f2.y = y2;
res.f2.s = sl;
res.f2.k = r_1;

    function [yy,a]=modelfun(k,x)
        a = m_1^k*(r_1-1)/(1-r_1*m_1^k);
        yy = r_1*(a+1).*(1-a./(a+x.^k));
    end

end