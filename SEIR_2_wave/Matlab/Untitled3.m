rng('default'); % For reproduceability.
pdTrue = GeneralizedGamma(1.4, 1.0, 1.6);
pdTrue = GeneralizedGamma(2.17525, 0.50125, 2.3655);
pdTrue = GeneralizedGamma(1/1.5625, 0.53125, 2.125/0.53125);
x=-1:0.01:10;
% pdf_x = shift/shape^rate*x.^(rate-1).*exp(-(x/shape).^shift)./gamma(rate/shift);
y = pdf(pdTrue,x);
% yr=real(y);
% close all;
% [~,ir]=min(yr);
% z0=[yr(ir:-1:1),yr(ir+1:end)];
% z0(ir-10:ir+35)=NaN;z=interp1(find(~isnan(z0)),z0(~isnan(z0)),1:length(x));
% z=z/sum(z);
plot(x,y);grid on;
% % hold on;
% % plot(x,z0);
% plot(x,z,'k');grid on;
n = 10000;
sample = pdTrue.drawSample(n);
pdEstimated = GeneralizedGamma();
pdEstimated.fitDist(sample)

