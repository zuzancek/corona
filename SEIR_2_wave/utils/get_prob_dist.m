function [opt_fit,kernel_fit,pdf_grid,cdf_grid,rv_grid] = get_prob_dist(N_rand,time_grid,pdf_base_grid,varargin)

ip = inputParser;
addParamValue(ip, 'do_plot', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'do_fitdist', true, @islogical);%#ok<*NVREPL>
addParamValue(ip, 'title', '', @ischar);%#ok<*NVREPL>
addParamValue(ip, 'dist_list', {'Kernel','Weibull','Gamma','Burr','GeneralizedExtremeValue','Lognormal','Loglogistic'}, @iscell);%#ok<*NVREPL>
parse(ip, varargin{:});
results = ip.Results;
do_plot = results.do_plot;
do_fitdist = results.do_fitdist;
dist_list = results.dist_list;
tit = results.title;

rnd_val = rand(N_rand,1);
min_val = min(rnd_val); max_val =  max(rnd_val);
rnd_val = (rnd_val-min_val)./(max_val-min_val);
% key variable
pdf_grid = time_grid.*pdf_base_grid; pdf_grid = pdf_grid./sum(pdf_grid);
cdf_grid = cumsum(pdf_grid);
idx = cell2mat(arrayfun(@(u) find(u<cdf_grid,1),rnd_val,'UniformOutput',false));
rv_grid = time_grid(idx);
kernel_fit = fitdist(rv_grid,'Kernel');
opt_fit = kernel_fit;

if do_plot
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
    dist = zeros(n);
    for i=1:n
        d{i}.type = dist_list{i};
        d{i}.obj = fitdist(rv_grid,dist_list{i});
        d{i}.cdf = cdf(d{i}.obj,time_grid);
        d{i}.pdf = pdf(d{i}.obj,time_grid);
        dist(i) = norm(pdf_grid-d{i}.pdf);
        d{i}.diff = dist(i);
        d{i}.mean = mean(d{i}.obj);
        d{i}.median = median(d{i}.obj);
        d{i}.std = std(d{i}.obj);
        d{i}.time_grid = time_grid;
        fprintf('%s error %2.4f\n',d{i}.type,d{i}.diff);
    end
    [~,idx] = sort(dist);
    idx = idx(2);
    opt_fit = d{idx};
    if do_plot
        figure('Name',strcat('Fitting Distribution (',tit,')'));
        col = [1 1 1];
        histogram(rv_grid',length(time_grid),'Normalization','probability','binwidth',1,...
            'FaceColor',0.65*col,'EdgeColor',0.5*col,'FaceAlpha',0.5,'EdgeAlpha',0.5); hold on;
        for i=1:n
            if i==idx
                plot(time_grid,d{i}.pdf/sum(d{i}.pdf),'linewidth',2,'color','g');
            elseif i==1                
                plot(time_grid,d{i}.pdf/sum(d{i}.pdf),'linewidth',1,'color','k');
            else
                plot(time_grid,d{i}.pdf/sum(d{i}.pdf),'linewidth',1);
            end
        end
        legend({'data',dist_list{:}}); %#ok<*CCAT>
        grid on;
    end
end

end

