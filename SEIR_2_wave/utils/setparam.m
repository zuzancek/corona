function [s] = setparam()

s = struct;
s.pop_size = 5443120;
s.sim_num = 10^6;
s.T = 200;

s.T_inf.mean = 2.9;
s.T_inf.std = 0.62;
s.T_inc.mean = 5.1;
s.Tinc.std = 0.62;
s.T_overlay = 1.5;
s.T_rem.mean = s.T_inf.mean+s.T_inc.mean-s.T_overlay;

g1.r0_scale = 10^2;
g1.mean = 0.98;
g1.alpha = 1;
s.g1 = g1;
g2.a0 = g1.mean+1/g1.r0_scale;
g2.a1 = 14;
g2.scale = 2.3;
g2.alpha = 1;
s.g2 = g2;
s.w_vec_default = 0.5+zeros(s.T,1);

s.obs_ratio = 1/16;

s.smooth_width = 3;
s.smooth_type = 3;
s.smooth_ends = 1;

end

