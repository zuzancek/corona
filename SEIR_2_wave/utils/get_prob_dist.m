function [pd_kernel,rv_grid] = get_prob_dist(N_rand,x_grid,cdf_grid,varargin)

ip = inputParser;
addParamValue(ip, 'do_plot', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'do_fitdist', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'dist_list', {'Kernel','Weibull','Poisson'}, @iscell);%#ok<*NVREPL>
parse(ip, varargin{:});
results = ip.Results;
do_plot = results.do_plot;
do_fitdist = results.do_fitdist;
dist_list = results.dist_list;

rnd_val = rand(N_rand,1);
% min_val = min(rnd_val); max_val =  max(rnd_val);
% rnd_val = (rnd_val-min_val)./(max_val-min_val);
% idx = cell2mat(arrayfun(@(u) (u==0).*1 +(u>0).*find(u<cdf_grid,1),rnd_val,'UniformOutput',false));
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

if do_fitdist
    n = length(dist_list);
    d = cell(n,1);
    for i=1:n
        d{i}.type = dist_list{i};
        d{i}.obj = fitdist(rv_grid,dist_list{i});
        d{i}.cdf = cdf(d{i}.obj,x_grid);
        d{i}.pdf = pdf(d{i}.obj,x_grid);
        d{i}.diff = norm(cdf_grid-d{i}.cdf);
        fprintf('%s error %2.4f\n',d{i}.type,d{i}.diff);
    end
    if do_plot
        figure;
        histogram(rv_grid',length(x_grid),'Normalization','probability','binwidth',1); hold on;
        for i=1:n
            plot(x_grid,d{i}.pdf/sum(d{i}.pdf),'linewidth',2);
        end
        legend({'data',dist_list{:}});
        grid on;
    end
end

end

