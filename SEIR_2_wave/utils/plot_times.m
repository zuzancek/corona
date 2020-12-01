function []=plot_times(s)
close all;
figure;
ay = 3;
T_lat = s.T_lat.mean;
T_pre = s.T_pre.mean;
T_test0 = s.T_test0;
T_inf_obs_0 = s.T_inf_obs0.mean;
T_inf_unobs = T_inf_obs_0+T_test0+T_pre;
ax = 1;
p1 = plot([0 T_lat/ax],[ay ay],'linewidth',1);hold on;
p2 = plot([T_lat/ax (T_lat+T_pre)/ax],[ay ay],'linewidth',1);hold on;
p3 = plot([(T_lat+T_pre)/ax (T_lat+T_pre+T_test0)/ax],[ay ay],'linewidth',1);hold on;
p4 = plot([(T_lat+T_pre+T_test0)/ax (T_lat+T_pre+T_test0+T_inf_obs_0)/ax],[ay ay],'linewidth',1);hold on;
p5 = plot([(T_lat)/ax (T_lat+T_pre+T_test0+T_inf_obs_0)/ax],[ay ay],'linewidth',1);hold on;
drawbrace([0 ay],[T_lat/ax ay],10,'Color',[1 1 1]);
text(.5*T_lat/ax,T_lat/ax+0.01,'T_{lat}','color',p1.Color);
drawbrace([0 ay],[T_lat/ax ay],10,'Color',p1.Color);
text((T_lat+0.5*T_pre)/ax,T_lat/ax+0.01,'T_{pre}','color',p2.Color);
drawbrace([T_lat/ax ay],[(T_lat+T_pre)/ax ay],10,'Color',p2.Color);
drawbrace([(T_lat+T_pre)/ax ay],[(T_lat+T_pre+T_test0)/ax ay],10,'Color',p3.Color);
text((T_lat+T_pre+0.5*T_test0)/ax,T_lat/ax+0.01,'T_{test}','color',p3.Color);
drawbrace([(T_lat+T_pre+T_test0)/ax ay],[(T_lat+T_pre+T_test0+T_inf_obs_0)/ax ay],10,'Color',p4.Color);
text((T_lat+T_pre+T_test0+0.5*T_inf_obs_0)/ax,T_lat/ax+0.01,'T_{inf}^{obs,0}','color',p4.Color);
drawbrace([(T_lat)/ax ay],[(T_lat+T_pre+T_test0+T_inf_obs_0)/ax ay],-10,'Color',p5.Color);
text((T_lat++0.5*(T_pre+T_test0+T_inf_obs_0))/ax,T_lat/ax-0.01,'T_{inf}^{unobs}','color',p5.Color);




