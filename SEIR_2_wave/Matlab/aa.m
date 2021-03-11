function []=simulate_model(s,init,dateFrom,dateTo)

T = dateTo-dateFrom+1;
burnin = s.firstData_offset;
firstData = -burnin+dateFrom;
rng(1000);

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
% Q_mat: transition times between states, size: Nt x N, different for o and y
create_Q_matrix();
% S_mat: current state for each agent, size Ns x N, different for o and y
create_S_matrix();
% P_mat: state-transition probabilities, size Ns x Ns, different for o and y
create_P_matrix();
% ST_mat: state-transition mapping, size Ns x Ns
% TS_mat: transition-state mapping, size Nt x 2
% U_mat: U-distr. random _vectors for transition probabilities in nodes with
% multiple out-paths, size 3 x N, different for o and y
create_aux_matrices();


%%
% reproduction number
Rt_vec = init.Rt(1:T_total);
betta_vec = Rt_vec./

% *** _vectors
V_vec = zeros(1,N); % vaccination
I_vec = zeros(1,N); % immunity (after infection)


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

    function []=get_E_inflow()
    end

    function []=create_Q_matrix()
        % transition time _matrix
        Q_mat_y(1,:) = 1;                                % SE
        Q_mat_y(2,:) = random(s.obj_lat,1,N);            % EU
        Q_mat_y(3,:) = random(s.obj_pre_test,1,N);       % UI
        Q_mat_y(4,:) = random(s.obj_pre_inf,1,N);        % URi
        Q_mat_y(5,:) = random(s.obj_ih_y,1,N);           % IH
        Q_mat_y(6,:) = random(s.obj_ir_y,1,N);           % IRi
        Q_mat_y(7,:) = random(s.obj_hd_y,1,N);           % HD
        Q_mat_y(8,:) = random(s.obj_hr_y,1,N);           % HRh
        Q_mat_y(9,:) = random(s.obj_im_i,1,N);           % RiS
        Q_mat_y(10,:) = random(s.obj_im_h,1,N);          % RhS
        Q_mat_o = Q_mat_y;
        Q_mat_o(5,:) = random(s.obj_ih_o,1,N);           % IH
        Q_mat_o(6,:) = random(s.obj_ir_o,1,N);           % IRi
        Q_mat_o(7,:) = random(s.obj_hd_o,1,N);           % HD
        Q_mat_o(8,:) = random(s.obj_hr_o,1,N);           % HRh
    end

    function []=create_S_matrix() %%  add init data here !!!!
        S_mat_y = zeros(Ns,N);
        S_mat_y(1,:) = 1;
        S_mat_o = zeros(Ns,N);
        S_mat_o(1,:) = 1;
        SR_mat_y = mod(find(S_mat_y==1)-1,Ns)+1;
        SR_mat_o = mod(find(S_mat_o==1)-1,Ns)+1;
    end

    function []=create_P_matrix()
        P_mat_y = zeros(Ns,Ns); % S E U I H D Ri Rh
        P_mat_y(1,2) = 1; P_mat_y(2,3) = 1; P_mat_y(7,1) = 1; P_mat_y(8,1) = 1; P_mat_o = P_mat_y;
        P_mat_y(3,4) = init.obs_y_ss; P_mat_y(3,7) = 1-init.obs_y_ss;
        P_mat_y(4,5) = init.eta_y_ss; P_mat_y(4,7) = 1-init.eta_y_ss;
        P_mat_y(5,6) = init.omega_y_ss; P_mat_y(5,8) = 1-init.omega_y_ss;
        P_mat_o(3,4) = init.obs_y_ss; P_mat_o(3,7) = 1-init.obs_y_ss;
        P_mat_o(4,5) = init.eta_y_ss; P_mat_o(4,7) = 1-init.eta_y_ss;
        P_mat_o(5,6) = init.omega_y_ss; P_mat_o(5,8) = 1-init.omega_y_ss;
        PR_mat_y = P_mat_y(3:6,:);
        PR_mat_o = P_mat_o(3:6,:);
    end

    function []=create_aux_matrices()
        U_mat_y = rand([3,N]);
        U_mat_o = rand([3,N]);  
        ST_mat = zeros(Ns,Ns);
        ST_mat(1,2)=1;        ST_mat(2,3)=2;        ST_mat(3,4)=3;         ST_mat(3,7)=4;         ST_mat(4,5)=5;
        ST_mat(4,7)=6;        ST_mat(5,6)=7;        ST_mat(5,8)=8;         ST_mat(7,1)=9;         ST_mat(8,1)=10;   
        TS_mat = zeros(Nt,2);
        TS_mat(1,1) = 1;       TS_mat(1,2) = 2;
        TS_mat(2,1) = 2;       TS_mat(2,2) = 3;
        TS_mat(3,1) = 3;       TS_mat(3,2) = 4;
        TS_mat(4,1) = 3;       TS_mat(4,2) = 7;
        TS_mat(5,1) = 4;       TS_mat(5,2) = 5;
        TS_mat(6,1) = 4;       TS_mat(6,2) = 7;
        TS_mat(7,1) = 5;       TS_mat(7,2) = 6;
        TS_mat(8,1) = 5;       TS_mat(8,2) = 8;
        TS_mat(9,1) = 7;       TS_mat(9,2) = 1;
        TS_mat(10,1) = 8;      TS_mat(10,2) = 1;                
        NT_mat = [3 5 7;4 6 8]; % p, 1-p
    end 
   
end