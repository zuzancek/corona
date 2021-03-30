
% t0=2.5;
% t1=8;
% ti0=1.5;
% ti1=3;
close all;
clear all;

x_max=15;
x=0:0.01:x_max;
m=4.975;s=1.175;
a=m*s*s;b=1/(s*s);
pd=makedist('Gamma','a',a,'b',b);
pdf_x=pdf(pd,x);
cdf_x=cdf(pd,x);
[~,idx_max]=max(pdf_x);
[~,idx_0] = min(abs(pdf_x(1:idx_max)-pdf_x(idx_max)/2));
[~,idx_2] = min(abs(pdf_x(idx_max:end)-pdf_x(idx_max)/2));
z = sum(pdf_x)*pdf_x.^.7/sum(pdf_x.^.7);
idx_2 = idx_2+idx_max;
t0 = x(idx_0);
ps =.45;
[~,idx_1] = min(abs(cdf_x-ps));
t1 = x(idx_1);
t2 = x(idx_2);
t_pre=t1-t0;
t_inf=t2-t0;
fprintf('T_pre: %2.4f\n',t_pre);
fprintf('T_inf: %2.4f\n',t_inf);
test = 2;
[~,t_idx] = min(abs(x(idx_1:end)-(t1+test)));
t_idx = t_idx+idx_1;
loss = cdf_x(t_idx);
fprintf('Loss (test=offset+2): %2.4f\n',loss);
cdf_x=cdf(pd,x);
subplot(2,1,1)
plot(x,pdf_x); xticks([0:x_max]);grid on;
hold on;
stem(x(idx_0),pdf_x(idx_0),'r');
stem(x(idx_2),pdf_x(idx_2),'r');
stem(x(idx_1),pdf_x(idx_1),'g');
stem(x(t_idx),pdf_x(t_idx),'k');
subplot(2,1,2)
plot(x,cdf_x); xticks([0:x_max]);grid on;
hold on;
stem(x(idx_0),cdf_x(idx_0),'r');
stem(x(idx_2),cdf_x(idx_2),'r');
stem(x(idx_1),cdf_x(idx_1),'g');
stem(x(t_idx),cdf_x(t_idx),'k');