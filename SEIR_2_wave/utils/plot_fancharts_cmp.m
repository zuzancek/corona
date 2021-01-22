function [] = plot_fancharts_cmp(q1_mat,q2_mat,s,disp_from,disp_to,varargin)

ip = inputParser;
addParamValue(ip, 'lab_idx', 1, @isnumeric); %#ok<*NVREPL>
addParamValue(ip, 'legend', {'A','B'}, @iscell);
addParamValue(ip, 'title','', @ischar);
addParamValue(ip, 'darkcenter',false, @islogical);
addParamValue(ip, 'offsetdate',dd(2020,9,1), @isnumeric);
parse(ip, varargin{:});
results = ip.Results;
lab_idx = results.lab_idx;
leg = results.legend;
tit = results.title;
offsetdate = results.offsetdate;

figure;
cm = lines(2);
idx_central1 = ceil(size(q1_mat,1)/2);
idx1 = (length(s.quant)-size(q1_mat,1))/2; quant1 = s.quant(idx1+1:end-idx1);
h1 = fanChart(1:size(q1_mat,2), q1_mat', q1_mat(idx_central1,:), quant1,...
    'alpha', .25, 'colormap', {'shadesOfColor',cm(1,:)},'midcolor',cm(1,:)/2,'darkcenter',true);
idx_central2 = ceil(size(q2_mat,1)/2);
idx2 = (length(s.quant)-size(q2_mat,1))/2; quant2 = s.quant(idx2+1:end-idx2);
h2 = fanChart(1:size(q2_mat,2), q2_mat', q2_mat(idx_central2,:),quant2,...
    'alpha', .25, 'colormap', {'shadesOfColor',cm(2,:)},'midcolor',cm(2,:)/2,'darkcenter',true);

legend([h1 h2],leg);

grid on;

dt = 1;
p0 = dt+(disp_from-offsetdate+1);
p1 = disp_to-offsetdate+1;
lab = {'marec','april','maj','jun','jul','august','september','oktober','november','december','januar'};
mt = [0 31 30 31 30 31 31 30 31 30 31];
mtt = cumsum(mt(lab_idx:end));
xticks(mtt)
xticklabels(lab(lab_idx:end));
xlim([p0 p1]);

title(tit);

end