function [] = plot_fanchart(q_mat,s,dt,disp_from,disp_to,t0)

f = figure;
fanChart(1:size(q_mat,2), q_mat', q_mat(s.quant_idx_central,:), s.quant,...
    'alpha', .75, 'colormap', {'shadesOfColor',s.color_graph});
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
mt = [0 31 30 31 30 31 31 30];
mtt = cumsum(mt);
xticks(mtt)
xticklabels({'marec','april','maj','jun','jul','august','september'})
xlim([p0 p1]);

end