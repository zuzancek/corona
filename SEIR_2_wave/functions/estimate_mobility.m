function []=estimate_mobility(mob,Rt,dateFrom,delay)

r_dateFrom = dateFrom;
r_dateTo = enddate(Rt);
r = resize(Rt,r_dateFrom:r_dateTo);
r_data = double(r);

m_dateFrom = dateFrom-delay;
m_dateTo = r_dateTo-delay;
m = resize(mob.orig/100,m_dateFrom:m_dateTo);
m_raw = double(resize(mob.raw/100,m_dateFrom:m_dateTo));
m_data = double(m);

[~,m_dmin] = min(m);
[~,r_dmin] = min(r);
assert(r_dmin-m_dmin==delay);

%
figure;
plot(m,'linewidth',1); hold on;
plot(r,'linewidth',1);
grid on;
legend({'mobility','Rt'});

%
% curve fitting: R(m) = (1/(a_0-m^3)-a_1)*c_2
% R(0) = 0; R(1) = 
m_idx = find(m_data<1,1); a_m = interp1([m_data(m_idx-1) m_data(m_idx)],[0 1],1);
r_1 = r_data(m_idx-1)*a_m+(1-a_m)*r_data(m_idx);
r_idx = find(r_data<1,1); a_r = interp1([r_data(r_idx-1) r_data(r_idx)],[0 1],1);
m_1 = m_data(r_idx-1)*a_r+(1-a_r)*m_data(r_idx);
a = m_1^3*(r_1-1)/(1-r_1*m_1^3);
figure;
x = 0.5:0.01:1.25;
y = r_1*(a+1).*(1-a./(a+x.^3));
plot(x,y,'linewidth',1);hold on;
scatter(m_data,r_data); hold on;
scatter(m_raw,r_data); hold on;
grid on;
xlabel('Mobility (%)');
ylabel('R_t');
legend({'fitted curve','smooth data','raw data'});

end