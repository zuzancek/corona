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
I0_y = init.I_y(firstData-shift_i+2:dateTo);
E0 = init.E(firstData-shift_i+2:dateTo);
E0_y = init.E_y(firstData-shift_i+2:dateTo);
H0 = init.H(firstData-shift_h+2:dateTo);
H0_y = init.H_y(firstData-shift_h+2:dateTo);
S0 = init.S(firstData-shift_h+2:dateTo);
S0_y = init.S_y(firstData-shift_h+2:dateTo);
D0 = init.D(firstData:dateTo);
D0_y = init.D_y(firstData:dateTo);
vac_state = init.vaccination.state;
vac_plan = init.vaccination.plan;

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
N_o = s.dep_ratio_65.*N;
N_y = N-N_o;
N_s = 8;
N_t = 10;
% *** matrices
%
% Q_mat: transition times between states, size: N_t x N, different for o and y
Q0_mat_y = zeros(N_t,N_y); Q0_mat_o = zeros(N_t,N_o);
Q_mat_y = zeros(N_t,N_y); Q_mat_o = zeros(N_t,N_o);
Q_mat_y_next = Q_mat_y; Q_mat_o_next = Q_mat_o;
create_Q_matrix();
%
% S_mat: current state for each agent, size N_s x N, different for o and y
% St_mat: time-evolution of states, aggregated, size T_total x N_s
S_mat_y = []; S_mat_o = []; SR_mat_y = []; SR_mat_o = []; St_mat_y=[]; St_mat_o=[]; St_mat=[];
S_mat_y_next = S_mat_y; S_mat_o_next = S_mat_o;
create_S_matrix();
%
% P_mat: state-transition probabilities, size N_s x N_s, different for o and y
P_mat_y = zeros(N_s,N_s); P_mat_o = zeros(N_s,N_s);
create_P_matrix();
%
% ST_mat: state-transition mapping, size N_s x N_s
% TS_mat: transition-state mapping, size N_t x 2
% U_mat: U-distr. random _vectors for transition probabilities in nodes with
% multiple out-paths, size 3 x N, different for o and y
U_mat_y = rand([3,N_y]);
U_mat_o = rand([3,N_o]);

%%
% reproduction number
Rt_vec = init.Rt(1:T_total);
T_inf_vec = random(s.obj_inf,1,N);
% transmission rate (per capita)
betta_vec = Rt_vec./T_inf_vec;
profile_U_y = []; profile_I_y = []; profile_U_o = []; profile_I_o = [];
SE_y_vec = zeros(1,N_y); RiS_y_vec = SE_y_vec; RhS_y_vec = SE_y_vec; EU_y_vec = SE_y_vec; IH_y_vec = SE_y_vec;
URi_y_vec = SE_y_vec; IRi_y_vec = SE_y_vec; HRh_y_vec = SE_y_vec; UI_y_vec = SE_y_vec; HD_y_vec = SE_y_vec;
SE_o_vec = zeros(1,N_o); RiS_o_vec = SE_o_vec; RhS_o_vec = SE_o_vec; UI_o_vec = SE_o_vec; IH_o_vec = SE_o_vec;
URi_o_vec = SE_o_vec; IRi_o_vec = SE_o_vec; HRh_o_vec = SE_o_vec; EU_o_vec = SE_o_vec; HD_o_vec = SE_o_vec;
get_inf_profile();

% *** _vectors
V_y_vec = zeros(1,N_y); V_o_vec = zeros(1,N_o); % vaccination
Ii_y_vec = zeros(1,N_y); Ii_o_vec = zeros(1,N_o); % immunity (after infection)
Ih_y_vec = zeros(1,N_y); Ih_o_vec = zeros(1,N_o); % immunity (after infection)


for t=1:T_total
    store_states_transitions();
    update_vaccination_immunity(); % todo: vymysliet vakcinaciu
    update_suspectible();
    update_exposed();
    update_unobserved();
    update_observed();
    update_hospitalized();
    update_dead();    
    update_removed();
    update_states_transitions(t); % todo
end

    function []=store_states_transitions()
        % transitions
        Q_mat_y_next = Q_mat_y; Q_mat_o_next = Q_mat_o;
        % states
        S_mat_y_next = S_mat_y; S_mat_o_next = S_mat_o;
    end

    function []=update_states_transitions(t)
        Q_mat_y = Q_mat_y_next; Q_mat_o = Q_mat_o_next;
        S_mat_y = S_mat_y_next; S_mat_o = S_mat_o_next;
        update_state_time(t);        
    end

%     function [v]=gen_random_pos(obj)
%         v = random(obj,1,N);
%         i0 = find(v<=0);
%         v(i0) = [];
%         cont=~isempty(i0);
%         while cont
%             v0 = random(obj,1,2*length(i0));
%             v = [v,v0];
%             v = v(v>0);
%             cont = length(v)<N;
%         end
%         v = v(1:N);
%     end

    function []=update_suspectible()
        % transitions: in 9,10 (RiS, RhS); out 1 (SE)
        % state: 1
        % calculate outflow of suspectible individuals (=inflow of exposed)
        calculate_E_inflow();
        % update transition matrices
        Q_mat_y_next(1,:) = Q_mat_y(1,:)-1;
        Q_mat_o_next(1,:) = Q_mat_o(1,:)-1;
        % newly suspectible time spent in this state is set to 1 and those
        % who are exposed have 0
        Q_mat_y_next(1,RiS_y_vec+RhS_y_vec) = T_total;
        Q_mat_o_next(1,RiS_o_vec+RhS_o_vec) = T_total;
        Q_mat_y_next(1,SE_y_vec) = 0;
        Q_mat_o_next(1,SE_o_vec) = 0;
        % update state
        S_mat_y_next(1,:) = S_mat_y_next(1,:)-SE_y_vec+RiS_y_vec+RhS_y_vec; % outflow of suspectible, new exposed
        S_mat_o_next(1,:) = S_mat_o_next(1,:)-SE_o_vec+RiS_o_vec+RhS_o_vec; % outflow of suspectible, new exposed
    end

    function []=update_exposed()
        % transitions: in 1 (SE), out 2 (EU)
        % state: 2
        % decrement time elapsed
        Q_mat_y_next(2,:) = Q_mat_y(2,:)-1;
        Q_mat_o_next(2,:) = Q_mat_o(2,:)-1;
        % determine agents leaving current state
        EU_y_vec = (Q_mat_y_next(2,:)==0);
        EU_o_vec = (Q_mat_o_next(2,:)==0);
        % newly exposed time spent in this state is set to their
        % corresponding def.values
        Q_mat_y_next(2,EU_y_vec) = Q0_mat_y(2,EU_y_vec);
        Q_mat_o_next(2,EU_o_vec) = Q0_mat_o(2,EU_o_vec);
        % update state
        S_mat_y_next(2,:) = S_mat_y_next(2,:)+SE_y_vec-EU_y_vec;
        S_mat_o_next(2,:) = S_mat_o_next(2,:)+SE_o_vec-EU_o_vec;
    end

    function []=update_unobserved()
        % transitions: in 2 (EU); out 3 (UI), 4 (URi)
        % state: 3
        % decrement time elapsed
        Q_mat_y_next(3,:) = Q_mat_y(3,:)-1;  Q_mat_y_next(4,:) = Q_mat_y(4,:)-1;
        Q_mat_o_next(3,:) = Q_mat_o(3,:)-1;  Q_mat_o_next(4,:) = Q_mat_o(4,:)-1;
        % determine agents leaving current state
        UI_y_vec = (Q_mat_y_next(3,:)==0);   URi_y_vec = (Q_mat_y_next(4,:)==0);
        UI_o_vec = (Q_mat_o_next(3,:)==0);   URi_o_vec = (Q_mat_o_next(4,:)==0);
        % newly infectious but yet unobserved (not tested/reported) time spent in this state is set to their
        % corresponding def.values
        % determine observed agents (transition 3), the remainder stays
        % uncovered (transition 4)
        [~,ui_new_y_idx] = maxk(U_mat_y(1,:).*EU_y_vec, ceil(P_mat_y(3,4)*sum(EU_y_vec))); 
        [~,ui_new_o_idx] = maxk(U_mat_o(1,:).*EU_o_vec, ceil(P_mat_o(3,4)*sum(EU_o_vec)));
        ur_new_y_idx = EU_y_vec-ui_new_y_idx;
        ur_new_o_idx = EU_o_vec-ui_new_o_idx;
        Q_mat_y_next(3,ui_new_y_idx) = Q0_mat_y(3,ui_new_y_idx);
        Q_mat_o_next(3,ui_new_o_idx) = Q0_mat_o(3,ui_new_o_idx);
        Q_mat_y_next(4,ur_new_y_idx) = Q0_mat_y(4,ur_new_y_idx);
        Q_mat_o_next(4,ur_new_o_idx) = Q0_mat_o(4,ur_new_o_idx);
        % update state
        S_mat_y_next(3,:) = S_mat_y_next(3,:)+EU_y_vec-UI_y_vec-URi_y_vec;
        S_mat_o_next(3,:) = S_mat_o_next(3,:)+EU_o_vec-UI_o_vec-URi_o_vec;
    end

    function []=update_observed()
        % transitions: in 3 (UI); out 5 (IH), 6 (IRi)
        % state: 4
        % decrement time elapsed
        Q_mat_y_next(5,:) = Q_mat_y(5,:)-1;  Q_mat_y_next(6,:) = Q_mat_y(6,:)-1;
        Q_mat_o_next(5,:) = Q_mat_o(5,:)-1;  Q_mat_o_next(6,:) = Q_mat_o(6,:)-1;
        % determine agents leaving current state
        IH_y_vec = (Q_mat_y_next(5,:)==0);   IRi_y_vec = (Q_mat_y_next(6,:)==0);
        IH_o_vec = (Q_mat_o_next(5,:)==0);   IRi_o_vec = (Q_mat_o_next(6,:)==0);
        % newly observed (new confirmed/reported cases) time spent in this state is set to their
        % corresponding def.values
        % determine agents who need hospitalization (transition 5), the
        % remainder is recovered at home (transition 6)
        [~,ih_new_y_idx] = maxk(U_mat_y(2,:).*UI_y_vec, ceil(P_mat_y(4,5)*sum(UI_y_vec))); 
        [~,ih_new_o_idx] = maxk(U_mat_o(2,:).*UI_o_vec, ceil(P_mat_o(4,5)*sum(UI_o_vec)));
        ir_new_y_idx = UI_y_vec-ih_new_y_idx;
        ir_new_o_idx = UI_o_vec-ih_new_o_idx;
        Q_mat_y_next(5,ih_new_y_idx) = Q0_mat_y(5,ih_new_y_idx);
        Q_mat_o_next(5,ih_new_o_idx) = Q0_mat_o(5,ih_new_o_idx);
        Q_mat_y_next(6,ir_new_y_idx) = Q0_mat_y(6,ir_new_y_idx);
        Q_mat_o_next(6,ir_new_o_idx) = Q0_mat_o(6,ir_new_o_idx);
        % update state
        S_mat_y_next(4,:) = S_mat_y_next(4,:)+UI_y_vec-IH_y_vec-IRi_y_vec;
        S_mat_o_next(4,:) = S_mat_o_next(4,:)+UI_o_vec-IH_o_vec-IRi_o_vec;
    end

    function []=update_hospitalized()
        % transitions: in 5 (IH); out 7 (HD), 8 (HRh)
        % state: 5
        % decrement time elapsed
        Q_mat_y_next(7,:) = Q_mat_y(7,:)-1;  Q_mat_y_next(8,:) = Q_mat_y(8,:)-1;
        Q_mat_o_next(7,:) = Q_mat_o(7,:)-1;  Q_mat_o_next(8,:) = Q_mat_o(8,:)-1;
        % determine agents leaving current state
        HD_y_vec = (Q_mat_y_next(7,:)==0);   HRh_y_vec = (Q_mat_y_next(8,:)==0);
        HD_o_vec = (Q_mat_o_next(7,:)==0);   HRh_o_vec = (Q_mat_o_next(8,:)==0);
        % newly hospitalized time spent in this state is set to their
        % corresponding def.values
        % determine agents who die (transition 7), the
        % remainder is recovered at hospital (transition 8)
        [~,hd_new_y_idx] = maxk(U_mat_y(3,:).*IH_y_vec, ceil(P_mat_y(5,6)*sum(IH_y_vec))); 
        [~,hd_new_o_idx] = maxk(U_mat_o(3,:).*IH_o_vec, ceil(P_mat_o(5,6)*sum(IH_o_vec)));
        hr_new_y_idx = IH_y_vec-hd_new_y_idx;
        hr_new_o_idx = IH_o_vec-hd_new_o_idx;
        Q_mat_y_next(7,hd_new_y_idx) = Q0_mat_y(7,hd_new_y_idx);
        Q_mat_o_next(7,hd_new_o_idx) = Q0_mat_o(7,hd_new_o_idx);
        Q_mat_y_next(8,hr_new_y_idx) = Q0_mat_y(8,hr_new_y_idx);
        Q_mat_o_next(8,hr_new_o_idx) = Q0_mat_o(8,hr_new_o_idx);
        % update state
        S_mat_y_next(5,:) = S_mat_y_next(5,:)+IH_y_vec-HD_y_vec-HRh_y_vec;
        S_mat_o_next(5,:) = S_mat_o_next(5,:)+IH_o_vec-HD_o_vec-HRh_o_vec;
    end

    function []=update_dead()
        % transitions: in 7 (IH); out - (final state)
        % state: 6
        % update state
        S_mat_y_next(6,:) = S_mat_y_next(6,:)+HD_y_vec;
        S_mat_o_next(6,:) = S_mat_o_next(6,:)+HD_o_vec;
    end

    function []=update_removed()
        % temporarily removed
        % transitions: in ; out 9 (RiS) and 10 (RhS)
        % states:  7 (Ri), 8 (Rh)
        % decrement time elapsed
        Q_mat_y_next(9,:) = Q_mat_y(9,:)-1;        Q_mat_y_next(10,:) = Q_mat_y(10,:)-1;
        Q_mat_o_next(9,:) = Q_mat_o(9,:)-1;        Q_mat_o_next(10,:) = Q_mat_o(10,:)-1;
        % newly removed time spent in this state is set to their
        % corresponding def.values
        Q_mat_y_next(9,URi_y_vec) = Q0_mat_y(9,URi_y_vec);
        Q_mat_y_next(10,HRh_y_vec) = Q0_mat_y(10,HRh_y_vec);
        Q_mat_y_next(9,IRi_y_vec) = Q0_mat_y(9,IRi_y_vec);
        Q_mat_o_next(9,URi_o_vec) = Q0_mat_o(9,URi_o_vec);
        Q_mat_o_next(10,HRh_o_vec) = Q0_mat_o(10,HRh_o_vec);
        Q_mat_o_next(9,IRi_o_vec) = Q0_mat_o(9,IRi_o_vec);
        % increase default removal time for those already temporarily
        % removed (better immunity)
        Q0_mat_y(9,URi_y_vec) = round(1.5*Q0_mat_y(9,URi_y_vec));
        Q0_mat_o(9,URi_o_vec) = round(1.5*Q0_mat_o(9,URi_o_vec));
        Q0_mat_y(9,iRi_y_vec) = round(1.5*Q0_mat_y(9,IRi_y_vec));
        Q0_mat_o(9,IRi_o_vec) = round(1.5*Q0_mat_o(9,IRi_o_vec));
        Q0_mat_o(10,HRh_y_vec) = round(1.5*Q0_mat_o(10,HRh_y_vec));
        Q0_mat_o(10,HRh_o_vec) = round(1.5*Q0_mat_o(10,HRh_o_vec));
        % newly suspectible (= outflow) - after loosing immunity
        RiS_y_vec = (Q_mat_y_next(9,:)==0);
        RhS_y_vec = (Q_mat_y_next(10,:)==0);
        RiS_o_vec = (Q_mat_o_next(9,:)==0);
        RhS_o_vec = (Q_mat_o_next(10,:)==0);
        % update state (+inflow -outflow)
        S_mat_y_next(7,:) = S_mat_y_next(7,:)+URi_y_vec+IRi_y_vec-RiS_y_vec;
        S_mat_y_next(8,:) = S_mat_y_next(8,:)+HRh_y_vec-RhS_y_vec;
        S_mat_o_next(7,:) = S_mat_o_next(7,:)+URi_o_vec+IRi_o_vec-RiS_o_vec;
        S_mat_o_next(8,:) = S_mat_o_next(8,:)+HRh_o_vec-RhS_o_vec;
    end

    function []=calculate_E_inflow()
        Fy = S_mat_y(1,:).*(s.psi_im_i.^Ii_y_vec)...
            .*(s.psi_im_h.^Ih_y_vec).*(s.psi_vac.^V_y_vec);
        Fo = s.alpha_s_o*S_mat_o(1,:).*(s.psi_im_i.^Ii_o_vec)...
            .*(s.psi_im_h.^Ih_o_vec).*(s.psi_vac.^V_o_vec);
        G = (sum(S_mat_y(3,:)+s.alpha_i_y.*S_mat_y(4,:)))./N_y ...
            + (sum(s.alpha_s_o.*S_mat_o(3,:)+s.alpha_i_o.*S_mat_o(4,:)))./N_o;
        Wy_vec = betta_vec.*Fy*G;
        Wo_vec = betta_vec.*Fo*G;
        ky = sum(Wy_vec);
        ko = sum(Wo_vec);
        [~,SE_y_idx] = maxk(Wy_vec,ky);
        [~,SE_o_idx] = maxk(Wo_vec,ko);
        SE_y_vec = 0*SE_y_vec; SE_y_vec(SE_y_idx) = 1;
        SE_o_vec = 0*SE_o_vec; SE_o_vec(SE_o_idx) = 1;
    end

    function []=update_vaccination_immunity()
        % immunity
        Ii_y_vec = Ii_y_vec.*(1+URi_y_vec+IRi_y_vec);
        Ih_y_vec = Ih_y_vec.*(1+HRh_y_vec);
        Ii_o_vec = Ii_o_vec.*(1+URi_o_vec+IRi_o_vec);
        Ih_o_vec = Ih_o_vec.*(1+HRh_o_vec);
        % vaccination
    end

    function []=create_Q_matrix()
        % transition time matrix (max per agent)
        Q0_mat_y(1,:) = T_total;                                   % SE
        Q0_mat_y(2,:) = round(random(s.obj_lat,1,N_y));            % EU
        Q0_mat_y(3,:) = round(random(s.obj_pre_test,1,N_y));       % UI
        Q0_mat_y(4,:) = round(random(s.obj_inf,1,N_y));            % URi
        Q0_mat_y(5,:) = round(random(s.obj_ih_y,1,N_y));           % IH
        Q0_mat_y(6,:) = round(random(s.obj_ir_y,1,N_y));           % IRi
        Q0_mat_y(7,:) = round(random(s.obj_hd_y,1,N_y));           % HD
        Q0_mat_y(8,:) = round(random(s.obj_hr_y,1,N_y));           % HRh
        Q0_mat_y(9,:) = round(random(s.obj_im_i,1,N_y));           % RiS
        Q0_mat_y(10,:) = round(random(s.obj_im_h,1,N_y));          % RhS
        Q0_mat_o(1,:) = T_total;                                   % SE
        Q0_mat_o(2,:) = round(random(s.obj_lat,1,N_o));            % EU
        Q0_mat_o(3,:) = round(random(s.obj_pre_test,1,N_o));       % UI
        Q0_mat_o(4,:) = round(random(s.obj_inf,1,N_o));            % URi
        Q0_mat_o(5,:) = round(random(s.obj_ih_o,1,N_o));           % IH
        Q0_mat_o(6,:) = round(random(s.obj_ir_o,1,N_o));           % IRi
        Q0_mat_o(7,:) = round(random(s.obj_hd_o,1,N_o));           % HD
        Q0_mat_o(8,:) = round(random(s.obj_hr_o,1,N_o));           % HRh
        Q0_mat_o(9,:) = round(random(s.obj_im_i,1,N_o));           % RiS
        Q0_mat_o(10,:) = round(random(s.obj_im_h,1,N_o));          % RhS
    end

    function []=create_S_matrix() %%  add init data here !!!!
        S_mat_y = zeros(N_s,N_y);
        S_mat_y(1,randperm(N_y,S0_y)) = 1;
        S_mat_y(2,randperm(N_y,E0_y)) = 1;
        S_mat_y(3,randperm(N_y,U0_y)) = 1;
        S_mat_y(4,randperm(N_y,I0_y)) = 1;
        S_mat_y(5,randperm(N_y,H0_y)) = 1;
        S_mat_y(6,randperm(N_y,D0_y)) = 1;
        S_mat_y(7,randperm(N_y,Ri0_y)) = 1;
        S_mat_y(8,randperm(N_y,Rh0_y)) = 1;
        S_mat_o = zeros(N_s,N_y);
        S_mat_o(1,randperm(N_y,S0-S0_y)) = 1;
        S_mat_o(2,randperm(N_y,E0-E0_y)) = 1;
        S_mat_o(3,randperm(N_y,U0-U0_y)) = 1;
        S_mat_o(4,randperm(N_y,I0-I0_y)) = 1;
        S_mat_o(5,randperm(N_y,H0-H0_y)) = 1;
        S_mat_o(6,randperm(N_y,D0-D0_y)) = 1;
        S_mat_o(7,randperm(N_y,Ri0-Ri0_y)) = 1;
        S_mat_o(8,randperm(N_y,Rh0-Rh0_y)) = 1;   
        update_state_time(1);
    end

    function []=update_state_time(t)
        SR_mat_y = mod(find(S_mat_y==1)-1,N_s)+1;
        SR_mat_o = mod(find(S_mat_o==1)-1,N_s)+1;
        St_mat_y = zeros(T_total,N_s);
        St_mat_y(t,:) = sum(SR_mat_y,2)';
        St_mat_o = zeros(T_total,N_s);
        St_mat_o(t,:) = sum(SR_mat_o,2)';
        St_mat = St_mat_y+St_mat_o;
    end

    function []=create_P_matrix()
        % S E U I H D Ri Rh
        P_mat_y(1,2) = 1; P_mat_y(2,3) = 1; P_mat_y(7,1) = 1; P_mat_y(8,1) = 1; P_mat_o = P_mat_y;
        P_mat_y(3,4) = init.xi_y; P_mat_y(3,7) = 1-init.xi_y;
        P_mat_y(4,5) = init.eta_y; P_mat_y(4,7) = 1-init.eta_y;
        P_mat_y(5,6) = init.omega_y; P_mat_y(5,8) = 1-init.omega_y;
        P_mat_o(3,4) = init.xi_o; P_mat_o(3,7) = 1-init.xi_o;
        P_mat_o(4,5) = init.eta_o; P_mat_o(4,7) = 1-init.eta_o;
        P_mat_o(5,6) = init.omega_o; P_mat_o(5,8) = 1-init.omega_o;
    end

    function []=get_inf_profile()
        mu = max(Q0_mat_y(4,:));
        profile_U_y = ones(ceil(mu),N_y);
        mi = max(Q0_mat_y(6,:));
        profile_I_y = ones(ceil(mi),N_y);
        mu = max(Q0_mat_o(4,:));
        profile_U_o = ones(ceil(mu),N_o);
        mi = max(Q0_mat_o(6,:));
        profile_I_o = ones(ceil(mi),N_o);
    end

end