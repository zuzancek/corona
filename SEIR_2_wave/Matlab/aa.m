function []=simulate_model(s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
burnin = s.firstData_offset;
firstData = -burnin+dateFrom;

T_total = dateTo-firstData+1;
shift_i = max([s.k_hosp,s.k_sick, s.k_ser]);
shift_h = max(s.k_death,s.k_rec);
shift = max(shift_i,shift_h);

I0 = init.I(firstData-shift_i+2:dateTo);
H0 = init.H(firstData-shift_h+2:dateTo);
S0 = init.S(firstData-shift_h+2:dateTo);
D0 = init.D(firstData:dateTo);

method_params = s.smoothing_method_params;
method_data = s.smoothing_method_data;

%% Equations
%
% .. V (vaccinated), R (reinfected), D (default, first infection, no vac.)
% . o (old), y (young)
%
% S_.._.(t) = S_.._.(t-1)+RS-SE_.._.(t);                            suspectible
% E_.._.(t) = E_.._.(t-1)+SE_.._.(t)-EU_.._.(t)                     exposed
% U_.._.(t) = U_.._.(t-1)+EU_.._.(t)-UR_.._.(t)-UM_.._.(t)          infectious, unobserved
% I_.._.(t) = I_.._.(t-1)+UI_.._.(t)-IR_.._.(t)-IH_.._.(t)          infectious, observed, home care
% H_.._.(t) = H_.._.(t-1)+IH_.._.(t)-HR_.._.(t)-HD_.._.(t)          hospital care
% D_.._.(t) = D_.._.(t-1)+HD_.._.(t)                                death
% R_.._.(t) = R_.._.(t-1)+UR_.._.(t)+OR_.._.(t)+HR_.._.(t)          recovery (temporal)
%

%%
N = s.pop_size;
Ns = 8;
Nt = 10;
% *** matrices
% Q: transition times between states, size: Nt x N, designed for o and y
create_Q_matrix();
% S: current state, size Ns x N, designed for o and y

% *** vectors
Vvec = zeros(1,N); % vaccination
Ivec = zeros(1,N); % immunity (after infection)


    function [v]=gen_random_pos(obj)
        v = random(obj,1,N);
        i0 = find(v<=0);
        v(i0) = [];
        cont=~isempty(i0);
        while cont
            v0 = random(obj,1,2*length(i0));
            v = [v,v0];
            v = v(v>0);
            cont = length(v)<N;
        end
        v = v(1:N);
    end

    function []=create_Q_matrix()
        % transition time matrix
        Qmat_y(1,:) = 1;                                % SE
        Qmat_y(2,:) = random(s.obj_lat,1,N);            % EU
        Qmat_y(3,:) = random(s.obj_pre_test,1,N);       % UI
        Qmat_y(4,:) = random(s.obj_pre_inf,1,N);        % URi
        Qmat_y(5,:) = random(s.obj_ih_y,1,N);           % IH
        Qmat_y(6,:) = random(s.obj_ir_y,1,N);           % IRi
        Qmat_y(7,:) = random(s.obj_hd_y,1,N);           % HD
        Qmat_y(8,:) = random(s.obj_hr_y,1,N);           % HRh
        Qmat_y(9,:) = random(s.obj_im_i,1,N);           % RiS
        Qmat_y(10,:) = random(s.obj_im_h,1,N);          % RhS
        Qmat_o = Qmat_y;
        Qmat_o(5,:) = random(s.obj_ih_o,1,N);           % IH
        Qmat_o(6,:) = random(s.obj_ir_o,1,N);           % IRi
        Qmat_o(7,:) = random(s.obj_hd_o,1,N);           % HD
        Qmat_o(8,:) = random(s.obj_hr_o,1,N);           % HRh
    end

    function []=create_S_matrix() %%  add init data here !!!!
        Smat_y = zeros(Ns,N);
        Smat_o = zeros(Ns,N);
        Pmat_y = zeros(Ns,Ns); % S E U I H D Ri Rh
        Pmat_y(1,2) = 1; Pmat_y(2,3) = 1; Pmat_y(7,1) = 1; Pmat_y(8,1) = 1; Pmat_o = Pmat_y;
        Pmat_y(3,4) = init.obs_y_ss; Pmat_y(3,7) = 1-init.obs_y_ss;
        Pmat_y(4,5) = init.eta_y_ss; Pmat_y(4,7) = 1-init.eta_y_ss;
        Pmat_y(5,6) = init.omega_y_ss; Pmat_y(5,8) = 1-init.omega_y_ss;
        Pmat_o(3,4) = init.obs_y_ss; Pmat_o(3,7) = 1-init.obs_y_ss;
        Pmat_o(4,5) = init.eta_y_ss; Pmat_o(4,7) = 1-init.eta_y_ss;
        Pmat_o(5,6) = init.omega_y_ss; Pmat_o(5,8) = 1-init.omega_y_ss;
    end

    

end
