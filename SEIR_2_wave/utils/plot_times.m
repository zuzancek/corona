function [] = plot_times(s)

close all;

drawArrow0 = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'Color','k',varargin{:});  

function [ h ] = drawarrow( x,y,~,~,props )
    %xlim(xlimits)
    %ylim(ylimits)
    h = annotation('arrow');
    set(h,'parent', gca, ...
        'position', [x(1),y(1),x(2)-x(1),y(2)-y(1)], ...
        'HeadLength', 10, 'HeadWidth', 10, 'HeadStyle', 'cback1', ...
        props{:} );
    drawArrow0(x,y,'linewidth',1);
end

figure('Name','COVID 19 characteristics');

ay = 3;
dy_plus = 0.2;
ymin = ay-0.1;
ymax = ay+dy_plus;
ep = 0.0025;
T_lat = s.T_lat.mean;
T_pre = s.T_pre.mean;
T_inc = T_lat+T_pre;
T_test0 = s.T_test0;
T_test = T_inc+T_test0;
T_inf_obs_0 = s.T_inf_obs0.mean;
SI = T_test+T_inf_obs_0;
xmax = ceil(SI);
% T_inf_unobs = T_inf_obs_0+T_test0+T_pre;
ax = 1;
p1 = plot([0 T_lat/ax],[ay ay],'linewidth',1);hold on;
p2 = plot([T_lat/ax (T_lat+T_pre)/ax],[ay ay],'linewidth',1);hold on;
p3 = plot([(T_lat+T_pre)/ax (T_lat+T_pre+T_test0)/ax],[ay ay],'linewidth',1);hold on;
p4 = plot([(T_lat+T_pre+T_test0)/ax (T_lat+T_pre+T_test0+T_inf_obs_0)/ax],[ay ay],'linewidth',1);hold on;
p5 = plot([(T_lat)/ax (T_lat+T_pre+T_test0+T_inf_obs_0)/ax],[ay ay]-ep,'linewidth',1);hold on;
plot([0 (T_lat+T_pre)/ax],[ay ay]+ep,'linewidth',1,'Color','m');hold on;

a = gca;
a.YLim = [ymin,ymax];
area([T_lat SI],[ymax ymax],'FaceColor',[130/255 230/255 0],'EdgeColor','k','FaceAlpha',0.25,'EdgeAlpha',0);

drawbrace([0 ay],[T_lat/ax ay],10,'Color',p1.Color);
text(.5*T_lat/ax,T_lat/ax+0.03,'T_{lat}','color',p1.Color);

drawbrace([T_lat/ax ay],[(T_lat+T_pre)/ax ay],10,'Color',p2.Color);
text((T_lat+0.5*T_pre)/ax,T_lat/ax+0.03,'T_{pre}','color',p2.Color);

drawbrace([(T_lat+T_pre)/ax ay],[(T_lat+T_pre+T_test0)/ax ay],10,'Color',p3.Color);
text((T_lat+T_pre+0.5*T_test0)/ax,T_lat/ax+0.03,'T_{test}','color',p3.Color);

drawbrace([(T_lat+T_pre+T_test0)/ax ay],[(T_lat+T_pre+T_test0+T_inf_obs_0)/ax ay],10,'Color',p4.Color);
text((T_lat+T_pre+T_test0+0.5*T_inf_obs_0)/ax,T_lat/ax+0.03,'T_{inf}^{obs,0}','color',p4.Color);

drawbrace([(T_lat)/ax ay],[(T_lat+T_pre+T_test0+T_inf_obs_0)/ax ay],-10,'Color',p5.Color);
text((T_lat+0.5*(T_pre+T_test0+T_inf_obs_0))/ax,T_lat/ax-0.03,'T_{inf}^{unobs,obs}','Color',p5.Color);

drawbrace([0 ay],[(T_lat+T_pre)/ax ay],-10,'Color','m');
text((0.5*(T_pre+T_lat))/ax,T_lat/ax-0.03,'T_{inc}','color','m');

drawarrow([T_inc T_inc],[ay+dy_plus/1.5 ay ],[0,xmax],[ymin,ymax],{'Color','k','LineWidth',1});
text(T_inc-0.5,ay+dy_plus/1.5+0.01,'Symptoms onset','color','k');

drawarrow([T_test T_test],[ay+dy_plus/2 ay ],[0,xmax],[ymin,ymax],{'Color','k','LineWidth',1});
text(T_test-0.125,ay+dy_plus/2+0.01,'Test','color','k');

ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';
ax.YMinorTick = 'off';
ax.YTick = [];
xlabel('days');

end




