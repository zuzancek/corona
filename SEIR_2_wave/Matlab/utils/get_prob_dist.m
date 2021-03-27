function [opt_fit,kernel_fit,pdf_grid,cdf_grid,rv_grid] = get_prob_dist(N_rand,time_grid,pdf_base_grid,varargin)

% input: time-dependent probability acounting for observation length 
% is calculated either using cumulated incidence 
% approach that accounts for competitive risk, or using Kaplan-Meier 
% (no need of competitive risk)

ip = inputParser;
addParamValue(ip, 'do_plot', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'do_ecdf', false, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'do_fitdist', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'censor', Inf, @isnumeric);%#ok<*NVREPL>
addParamValue(ip, 'title', '', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'tol',.05, @isnumeric);%#ok<*NVREPL>
addParamValue(ip, 'plot_num',5, @isnumeric);%#ok<*NVREPL>
addParamValue(ip, 'dist_list', ...
    {'Kernel','Weibull','Gamma','Burr','GeneralizedExtremeValue','InverseGaussian',...
    'Lognormal','Loglogistic','BirnbaumSaunders','Exponential','HalfNormal',...
    'Logistic','Nakagami','Rayleigh','GeneralizedPareto'}, @iscell);%#ok<*NVREPL>

parse(ip, varargin{:});
results = ip.Results;
do_plot = results.do_plot;
do_fitdist = results.do_fitdist;
dist_list = results.dist_list;
tit = results.title;
tol = results.tol;
N = results.plot_num;
do_ecdf = results.do_ecdf;
LargeNum = 1000;

rnd_val = rand(N_rand,1);
min_val = min(rnd_val); max_val =  max(rnd_val);
rnd_val = (rnd_val-min_val)./(max_val-min_val);

% key variable
pdf_grid = time_grid.*pdf_base_grid; 
pdf_grid = pdf_grid./sum(pdf_grid);
cdf_grid = cumsum(pdf_grid);
idx = cell2mat(arrayfun(@(u) find(u<cdf_grid,1),rnd_val,'UniformOutput',false));
rv_grid = time_grid(idx);
kernel_fit = fitdist(rv_grid,'Kernel');
opt_fit = kernel_fit;
f=ecdf(rv_grid); 
idx_last = find(f>1-tol,1);

if do_ecdf
    figure('Name',strcat('Empirical Distribution (',tit,')'));
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
    dist = zeros(n,1);
    for i=1:n
        d{i}.type = dist_list{i};
        try
            d{i}.obj = fitdist(rv_grid,dist_list{i});
            d{i}.cdf = cdf(d{i}.obj,time_grid);
            d{i}.pdf = pdf(d{i}.obj,time_grid); d{i}.pdf = d{i}.pdf./sum(d{i}.pdf);
            dist(i) = norm(pdf_grid(1:idx_last)-d{i}.pdf(1:idx_last));
            d{i}.diff = dist(i);
            d{i}.mean = mean(d{i}.obj);
            d{i}.median = median(d{i}.obj);
            d{i}.std = std(d{i}.obj);
            d{i}.time_grid = time_grid;
            fprintf('%s error %2.4f\n',d{i}.type,d{i}.diff);
        catch err %#ok<NASGU>
            d{i}.diff = LargeNum;
            dist(i) = LargeNum;
        end
    end
    [~,idx] = sort(dist);
    idx_ker = find(strcmp(dist_list,'Kernel'));
    idx_opt = idx(1);
    if idx_ker == idx_opt
        idx_opt = idx(2);
    end
    opt_fit = d{idx_opt};
    opt_fit.epdf = d{idx_ker}.pdf;
    opt_fit.emean = d{idx_ker}.mean;
    opt_fit.ecdf = d{idx_ker}.cdf;
    opt_fit.eobj = d{idx_ker}.obj;
    fprintf('\n********** Error minimizing distribution: %s\n',opt_fit.type);
    if do_plot
        figure('Name',strcat('Fitting Distribution (',tit,')'));
        col = [1 1 1];
        histogram(rv_grid',length(time_grid),'Normalization','probability','binwidth',1,...
            'FaceColor',0.65*col,'EdgeColor',0.5*col,'FaceAlpha',0.5,'EdgeAlpha',0.5); hold on;
        for i=1:N
            if i==find(idx==idx_opt) 
                plot(time_grid,d{idx_opt}.pdf/sum(d{idx_opt}.pdf),'linewidth',2,'color','m');
            elseif i==find(idx==idx_ker)               
                plot(time_grid,d{idx_ker}.pdf/sum(d{idx_ker}.pdf),'-.','linewidth',1,'color','k');
            else
                plot(time_grid,d{idx(i)}.pdf/sum(d{idx(i)}.pdf),'linewidth',1);
            end
        end
        legend({'data',dist_list{idx(1:N)}}); %#ok<*CCAT>
        grid on;
    end
end

end

