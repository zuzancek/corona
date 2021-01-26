function []=get_prob_die()

p_d_n = 2.71/100;   p_d_c = 4.92/100;   p_d_v = 26.21/100;
p_no = 31.98/100;   p_cn_o = 26/100;    p_co = p_no*p_cn_o; p_vc_o = 0.26;
p_ny = 2.32/100;    p_cn_y = 20/100;    p_cy = p_ny*p_cn_y; p_vc_y = 0.2;
p_d_y = 5.15/100;   p_d_y_cv = 15.75/100;
p_d_o = 37.16/100;  p_d_o_cv = 45.01/100;

    function [d] = solve_prob(x)
        d = zeros(6,1);
        d(1) = x(1)*p_no+x(2)*p_ny-p_d_n;
        d(2) = x(3)*p_co+x(4)*p_cy-p_d_c;
        d(3) = (x(5)*p_vc_o*p_co+x(6)*p_cy*p_vc_y)-p_d_v;
        d(4) = x(1)+x(3)*p_cn_o+x(5)*p_vc_o*p_cn_o-p_d_o;
        d(5) = x(2)+x(4)*p_cn_y+x(6)*p_vc_y*p_cn_y-p_d_y;
        d(6) = x(3)+x(5)*p_vc_o-p_d_o_cv;
        % d(7) = x(4)+x(6)*x(7)-p_d_y_cv;
    end

x0 = [0.03,0.005,0.065,0.02,0.7,0.1];
[a,b] = fsolve(@solve_prob,x0);

end