function [Rt,q_mat,res,Rt_last,Rt_dist,Rt_rnd] = estimate_Rt_SEIR(inputs,s,do_quant,do_weight,do_dist)

%% Model
% S(t+1) = S(t)-Z(t)
% E(t+1) = E(t)+Z(t)-E(t)/T_lat
% Iu(t+1) = Iu(t)+E(t)/T_lat-rho(t)*Iu(t)/T_test-(1-rho(t))*Iu(t)/T_inf
% Io(t+1) = Io(t)+rho(t)*Iu(t)/T_test-lambda(t)*Io(t)/T_hosp-(1-lambda(t))*Io(t)/T_sick
% H(t+1) = H(t)+lambda(t)*Io(t)/T_hosp-omega(t)*H(t)/T_death-(1-omega(t))*H(t)/T_rec
% D(t+1) = D(t)+omega(t)*H(t)/T_death
%
% Z(t) = Rt(t)/T_inf*S(t)/N*(Iu(t)+alpha_o*Io(t)+alpha_h*H(t))
% X(t) = rho(t)*Iu(t)/T_test

%% initialization
N = s.sim_num;
pop_size = s.pop_size;
X = inputs.Z;
T = length(X);
init = inputs.init;
params = inputs.params;
I0 = init.I;
H0 = init.H;
D0 = init.H;

try    
    rho = params.obs_ratio;
    assert(length(rho)>=T);
catch err %#ok<NASGU>
    rho = s.obs_ratio+0*X;
end
try    
    omega = params.old_ratio;
    assert(length(omega)>=T);
catch err %#ok<NASGU>
    omega = s.old_share+0*X;
end
try 
    sigma = params.asymp_ratio;
    assert(length(sigma)>=X);
catch err %#ok<NASGU>
    sigma = 1-s.symp_ratio_obs+0*Z;
end

T_test = s.T_test; % zmenit na T_pre (rv) + 1-2 dni postupne
T_inf = s.T_inf; T_inf_vec = get_rv(T_inf);
T_lat = s.T_lat; T_lat_vec = get_rv(T_lat);
T_sick_y.mean = s.T_sick_y.mean-T_test.mean; T_sick_y_vec = get_rv(T_sick_y);
T_sick_o.mean = s.T_sick_o.mean-T_test.mean; T_sick_o_vec = get_rv(T_sick_o);
alpha_ihy = s.eta_y/s.T_hosp_y;         alpha_iho = s.eta_o/s.T_hosp_o; 
alpha_iry = (1-s.eta_y)./T_sick_y_vec;  alpha_iro = (1-s.eta_o)./T_sick_o_vec;
theta_ih = omega.*alpha_iho+(1-omega).*alpha_ihy;
theta_ir = omega'.*alpha_iro+(1-omega').*alpha_iry;
varrho = omega.*alpha_iho./theta_ih;
alpha_hdy = s.omega_y/s.T_death_y;    alpha_hdo = s.omega_o/s.T_death_o;
alpha_hry = (1-s.omega_y)./p.T_rec_y; alpha_hro = (1-s.omega_o)./p.T_rec_o; 
theta_hd = varrho.*alpha_hdo+(1-varrho).*alpha_hdy;
theta_hr = varrho.*alpha_hro+(1-varrho).*alpha_hry;
theta_ui = rho./T_test.mean;
theta_ur = (1-rho')./T_inf_vec;

xx = 0:0.001:10;
yy = cdf('Gamma',xx,T_inf.mean*(T_inf.std^2),1/(T_inf.std^2));
zeta_o = 1-yy((find(xx>=T_test.mean,1)));                  alpha_o = 0.5*zeta_o;
zeta_h = 1-yy((find(xx>=T_test.mean+T_hosp.mean,1)));      alpha_h = 0.25*zeta_h;

% define arrays
S = zeros(N,T); E = S; Io = S; Iu = S; Rt = S; H = S; D = S;
S(:,1) = pop_size-I0; Io(:,1) = I0.*rho(1); Iu(:,1) = I0.*(1-rho(1));
H(:,1) = H0; D(:,1) = D0;

%%
idx = ones(N,1);

for t=1:T-2
    Io(:,t+1) = Io(:,t).*(1-theta_ih(t)-theta_ir(:,t))+X(t);
    H(:,t+1) = H(:,t).*(1-theta_hd(t)-theta_hr(t))+theta_ih(t).*Io(:,t);
    D(:,t+1) = D(:,t)+theta_hd(t).*H(:,t);
    Iu(:,t) = X(t)./theta_ui(t);
    Iu(:,t+1) = X(t+1)./theta_ui(t+1);
    Iu(:,t+2) = X(t+2)./theta_ui(t+2);
    E(:,t) = (Iu(:,t+1)-Iu(:,t).*(1-theta_ui(t)-theta_ur(t))).*T_lat_vec;
    E(:,t+1) = (Iu(:,t+2)-Iu(:,t+1).*(1-theta_ui(t+1)-theta_ur(t+1))).*T_lat_vec;
    Z = E(:,t+1)-E(:,t).*(1-1./T_lat_vec);
    Rt(:,t) = pop_size.*Z./S(:,t).*T_inf_vec./(Iu(:,t)+alpha_o.*Io(:,t)+alpha_h.*H(:,t));
    S(:,t+1) = S(:,t)-Z;        
    idx = idx & Io(:,t+1)>=0 & Iu(:,t+1)>=0 & H(:,t+1)>=0 & E(:,t+1)>=0;
end
idx = find(idx>0);
S = S(idx,:); E = E(idx,:); Iu = Iu(idx,:); Io = Io(idx,:); H = H(idx,:); D = D(idx,:); Rt = Rt(idx,:);
I = Iu+Io;
Ia = Io.*sigma; Is = Io-Ia;

%% store results
% means
for t = 1:T
    res.mean.Rt_mean(t) = mean(Rt(:,t));
    res.mean.I_mean(t) = mean(I(:,t));
    res.mean.Iu_mean(t) = mean(Iu(:,t));
    res.mean.Io_mean(t) = mean(Io(:,t));
    res.mean.Ia_mean(t) = mean(Ia(:,t));
    res.mean.Is_mean(t) = mean(Is(:,t));
    res.mean.S_mena(t) = mean(S(:,t));
    res.mean.H_mean(t) = mean(H(:,t));
    res.mean.D_mean(t) = mean(D(:,t));
    res.mean.E_mean(t) = mean(E(:,t));
end

%%

if do_weight
    weights = s.pweight;
    last_num = length(weights);
    nn = size(Rt_vec(:,T),1);
    weights_mat = repmat(weights,nn,1);
    Rt_last = Rt_vec(:,T-last_num+1:T);
    Rt_last = sum(Rt_last.*weights_mat,2);
else
    Rt_last = Rt_vec(:,T);
end

if do_quant
    q_vec = s.quant;
    M = length(q_vec);
    q_mat = zeros(M,T);
    for j = 1:M
        q_mat(j,:) = quantile(Rt_vec,q_vec(j),1);
    end
else
    q_mat = [];
end

if do_dist
    Rt_dist = cell(T,1);
    Rt_rnd = zeros(N,T);
    tol = 0.005;
    % create time-varying estimator of Rt
    for t = 1:T
        r_min = min(Rt_vec(:,t));
        r_max = max(Rt_vec(:,t));
        num_pts = s.min_pts+(min(s.max_dif,ceil(max(s.min_dif,r_max-r_min))-s.min_dif)/(s.max_dif-s.min_dif)*(s.max_pts-s.min_pts));
        pts = linspace(r_min,r_max,num_pts);
        % [cdf_x,x] = ksdensity(Rt_vec(:,t),pts,'Support','positive','BoundaryCorrection','reflection','Function','cdf');
        [cdf_x,x] = ksdensity(Rt_vec(:,t),pts,'Function','cdf');
        idx = find(cdf_x(2:end)-cdf_x(1:end-1)<=0,1);
        if ~isempty(idx)
            x = x(1:idx); cdf_x = cdf_x(1:idx);
        end
        if cdf_x(1)>=tol/2
            x = [x(1)-(x(2)-x(1)),x(1)-0.5*(x(2)-x(1)),x];
            cdf_x = [0,(cdf_x(2)-cdf_x(1))/(cdf_x(3)-cdf_x(1))*cdf_x(1),cdf_x]; %#ok<*AGROW>
        end
        Rt_dist{t} = {x,cdf_x};
        cutoff = max([1e-5,cdf_x(1),1-cdf_x(end)]);
        u_r = random('uniform',0+cutoff,1-cutoff,N,1);
        try
            x_r = interp1(cdf_x,x,u_r,'pchip','extrap');
        catch err %#ok<NASGU>
            disp(t);
        end
        Rt_rnd(:,t) = reshape(x_r,[],1);
    end
else
    Rt_dist = [];
    Rt_rnd = [];
end

    function [x] = get_rv(y)
        shape0 = y.mean.*(y.std)^2; scale0 = 1./(y.std)^2;
        L = length(shape0);
        shape0_vec = repmat(shape0,N,1);
        scale0_vec = scale0*ones(N,L);
        x = gamrnd(shape0_vec,scale0_vec);
    end

end