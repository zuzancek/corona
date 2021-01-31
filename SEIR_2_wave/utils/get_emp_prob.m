function [pd_kernel,rv_grid] = get_emp_prob(N_rand,x_grid,cdf_grid,do_plot)

rnd_val = rand(N_rand,1);
idx = cell2mat(arrayfun(@(u) find(u<cdf_grid,1),rnd_val,'UniformOutput',false));
rv_grid = x_grid(idx);
pd_kernel = fitdist(rv_grid,'Kernel');

if do_plot
    figure;
    [f,x]=ecdf(rv_grid); x(1)=0; 
    subplot(3,1,1); plot(x,f);grid on;title('CDF(x)');
    [f,x]=ecdf(rv_grid,'function','cumulative hazard'); x(1)=0;
    subplot(3,1,2); plot(x,f);grid on;title('H(x)');
    [f,x]=ecdf(rv_grid,'function','survivor'); x(1)=0;
    subplot(3,1,3); plot(x,f);grid on;title('S(x)');
end

end

