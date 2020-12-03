function [] = plot_fanchart(q_mat,s,dt,disp_from,disp_to,t0,fname,maketit,varargin)

f = figure('Name',fname);
fanChart(1:size(q_mat,2), q_mat', q_mat(s.quant_idx_central,:), s.quant,...
    'alpha', .75, 'colormap', {'shadesOfColor',s.color_graph});
if isempty(varargin)
    lab_idx = 1;
else
    lab_idx = varargin{1};
    if length(varargin)>1
        ydata = varargin{2};
        hold on;
        pp = plot(ydata.data,'color',ydata.col,'linewidth',2);
        legend(pp,ydata.leg);
    end
end

grid on;
set(gca,'color',s.color_bkg);
ax = gca;
ax.GridColor = s.color_grid;
ax.FontName = 'TimesNewRoman';
ax.FontWeight = 'bold';
ax.XAxis.Color = s.color_grid;
ax.YAxis.Color = s.color_grid;
f.Color = s.color_bkg;
p0 = dt+(disp_from-t0+1);
p1 = disp_to-t0+1;
lab = {'marec','april','maj','jun','jul','august','september','oktober','november','december'};
mt = [0 31 30 31 30 31 31 30 31 30 31];
mtt = cumsum(mt(lab_idx:end));
xticks(mtt)
xticklabels(lab(lab_idx:end));
xlim([p0 p1]);
if maketit
    title(fname,'Color',s.color_grid);
end

end