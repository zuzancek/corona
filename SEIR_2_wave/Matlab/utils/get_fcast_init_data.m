function [q]=get_fcast_init_data(p,T)

q.S_o = p.S_o(T);    q.S_y = p.S_y(T);
q.E_o = p.E_o(T);    q.E_y = p.E_y(T);
q.O_o = p.S_o(T);    q.O_y = p.O_y(T);
q.U_o = p.S_o(T);    q.U_y = p.U_y(T);
q.Rt_avg = mean(resize(p.Rt,T-7:T));
q.sigma_o_avg = mean(resize(p.sigma_o,T-7:T));
q.sigma_y_avg = mean(resize(p.sigma_y,T-7:T));

end