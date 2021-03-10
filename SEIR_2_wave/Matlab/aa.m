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
% O_.._.(t) = O_.._.(t-1)+UO_.._.(t)-OR_.._.(t)-OH_.._.(t)          infectious, observed, home care
% H_.._.(t) = H_.._.(t-1)+OH_.._.(t)-HR_.._.(t)-HD_.._.(t)          hospital care
% D_.._.(t) = D_.._.(t-1)+HD_.._.(t)                                death
% R_.._.(t) = R_.._.(t-1)+UR_.._.(t)+OR_.._.(t)+HR_.._.(t)          recovery (temporal)
%

%%
N = s.pop_size;
Ns = 8;
Nt = 10;
% *** matrices
% Q: 
% transition times between states, size: Nt x N
% ordering: SE,EU,UO,URd,OH,ORd,HD,HRh,RdS,RhS
Q = zeros(Nt,N);
Q(1,:) = 1;


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
        Q = zeros(Nt,N);
        
    end

end
