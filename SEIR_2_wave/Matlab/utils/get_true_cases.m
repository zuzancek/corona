function[] = get_true_cases()

s=setparam();
close all;

x_max_mult = 3.5;
x=0:0.01:x_max_mult*s.T_inf.mean;
pd_obj = s.obj_inf;
pdf_x = pdf(pd_obj,x);
cdf_x = cdf(pd_obj,x);

test_day = s.T_test.mean;
onset_day = s.T_pre.mean;
[~,idx_test_d2] = min(abs(x-(test_day+onset_day)));
[~,idx_test_d3] = min(abs(x-(test_day+1+onset_day)));
[~,idx_onset] = min(abs(x-(onset_day)));

figure('Name','Infectious period');
subplot(2,1,1)
plot(x,100*pdf_x,'linewidth',1); hold on;
p0=stem(x(idx_onset),100*pdf_x(idx_onset));
p1=stem(x(idx_test_d2),100*pdf_x(idx_test_d2),'g');
p2=stem(x(idx_test_d3),100*pdf_x(idx_test_d3),'c');
grid on;
ylabel('%');xlabel('days');
legend([p0 p1 p2],{'symptoms onset' 'testing (D=2)', 'testing (D=3)' });
title('Infectious period: PDF');
xticks(0:ceil(x_max_mult*s.T_inf.mean));
subplot(2,1,2)
plot(x,100*cdf_x,'linewidth',1);hold on;
p0=stem(x(idx_onset),100*cdf_x(idx_onset));
p1=stem(x(idx_test_d2),100*cdf_x(idx_test_d2),'g');
p2=stem(x(idx_test_d3),100*cdf_x(idx_test_d3),'c');
grid on;
ylabel('%');xlabel('days');
legend([p0 p1 p2],{'symptoms onset' 'testing (D=2)', 'testing (D=3)' });
grid on;
xticks(0:ceil(x_max_mult*s.T_inf.mean));
ylabel('%');xlabel('days');
title('Infectious period: CDF');

pow=.7;
pdf_i = sum(pdf_x)*pdf_x.^pow/sum(pdf_x.^pow);
cdf_i = cumsum(pdf_i)/sum(pdf_x);
cps=.45;
[~,idx_0] = min(abs(cdf_i-cps));
figure;
subplot(2,1,1);
plot(x,pdf_i);hold on;
stem(x(idx_0),pdf_i(idx_0));
grid on;
xticks(0:ceil(x_max_mult*s.T_inf.mean));
subplot(2,1,2);
plot(x,cdf_i);hold on;
stem(x(idx_0),cdf_i(idx_0));
grid on;
xticks(0:ceil(x_max_mult*s.T_inf.mean));
