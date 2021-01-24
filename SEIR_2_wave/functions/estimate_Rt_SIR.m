function [Rt,q_mat,res,Rt_last,Rt_dist,Rt_rnd] = estimate_Rt_SIR(inputs,s,do_quant,do_weight,do_dist)

%% Model
% S(t+1) = S(t)-Z(t)
% Iu(t+1) = Iu(t)+Z(t)-rho(t)*Iu(t)/T_test-(1-rho(t))*Iu(t)/T_inf
% Io(t+1) = Io(t)+rho(t)*Iu(t)/T_test-lambda(t)*Io(t)/T_hosp-(1-lambda(t))*Io(t)/T_sick
% H(t+1) = H(t)+lambda(t)*Io(t)/T_hosp-omega(t)*H(t)/T_death-(1-omega(t))*H(t)/T_rec
% D(t+1) = D(t)+omega(t)*H(t)/T_death
%%

% params.obs_ratio = cases_data.obs_ratio;           
% params.old_ratio = cases_data.old_ratio_smooth;
% params.death_ratio = deaths_data.old_ratio_smooth;
% params.asymp_ratio = cases_data.asymp_ratio_smooth;

% initialization
z_obs = inputs.z;
T = length(z_obs);
try
    obs_ratio = inputs.obs_ratio;
    obs_ratio = obs_ratio(1:T);
catch err %#ok<NASGU>
    obs_ratio = s.obs_ratio+0*inputs.z;
end
try
    sigma = inputs.asymp_ratio;
    assert(length(sigma)>=length(inputs.z));
    sigma = sigma(1:T);
catch err %#ok<NASGU>
    if isempty(inputs.asymp_ratio)
        sigma = (1-s.symp_ratio_obs)+0*inputs.z;
    elseif length(inputs.asymp_ratio)<length(inputs.z)
        nn = length(inputs.z);
        sigma(end:end+nn-length(sigma)) = sigma(end);
    end
end
try
    T_hosp = inputs.T_hosp;
    assert(length(T_hosp)>=length(inputs.z));
catch err %#ok<NASGU>
    T_hosp = s.T_hosp+0*inputs.z;
end
try
    rho = double(inputs.old_ratio);    
    assert(length(rho)>=length(inputs.z));
    rho = rho(1:T);
catch err %#ok<NASGU>
    rho = s.old_share+0*inputs.z;
end
try
    varsigma = double(inputs.death_ratio);    
    assert(length(varsigma)>=length(inputs.z));
catch err %#ok<NASGU>
    varsigma = s.old_death_ratio+0*inputs.z;
end
z = z_obs./obs_ratio(1:length(z_obs));
z_unobs = z-z_obs;
z_obs_o = z_obs.*rho;
z_obs_y = z_obs-z_obs_o;

I0 = inputs.I0/obs_ratio(1);
N = s.sim_num;
pop_size = s.pop_size;
s.SI.mean = 6.5;
T_si_vec = get_rv(s.SI);
% eta_hr = s.eta_hr;
% lambda = s.lambda; % 0.0579;
% T_death = s.T_death.mean; % 7.17; 
% T_rec = s.T_rec;
alpha_ihy = s.eta_y/s.T_hosp_y;        
alpha_iho = s.eta_o/s.T_hosp_o;  
alpha_iry = (1-s.eta_y)./(s.T_sick_y-s.T_test.mean); 
alpha_iro = (1-s.eta_o)./(s.T_sick_o-s.T_test.mean);
alpha_hdy = s.omega_y/s.T_death_y;
alpha_hdo = s.omega_o/s.T_death_o;
alpha_hry = (1-s.omega_y)./s.T_rec_y; 
alpha_hro = (1-s.omega_o)./s.T_rec_o; 
% T_hosp = s.T_hosp.mean;
alpha = 1.175; %((s.T_inf_obs.mean-s.T_inf_obs0.mean)+s.T_inf_obs0.mean/s.case_isolation_effect)/s.T_inf_unobs.mean;
% set initial values
S_vec = zeros(N,T); S_vec(:,1) = pop_size-I0;
I_vec = zeros(N,T); I_vec(:,1) = I0;
I_obs_vec = zeros(N,T); I_obs_vec(:,1) = I0*obs_ratio(1);
I_unobs_vec = zeros(N,T); I_unobs_vec(:,1) = I0-I_obs_vec(:,1);
I_obs_o_vec = zeros(N,T); I_obs_o_vec(:,1) = I0*s.old_share;
I_obs_y_vec = zeros(N,T); I_obs_y_vec(:,1) = I0-I_obs_o_vec(:,1);
% I_asympt_vec = zeros(N,T); I_asympt_vec(:,1) = I_obs_vec(:,1).*sigma(1);
% I_sympt_vec = zeros(N,T); I_sympt_vec(:,1) = I_obs_vec(:,1).*(1-sigma(1));
H_vec = zeros(N,T); H_o_vec = zeros(N,T); H_y_vec = zeros(N,T); 
F_vec = zeros(N,T); 
D_vec = zeros(N,T); D_o_vec = zeros(N,T); D_y_vec = zeros(N,T); 
if isfield(inputs, 'H0')
    H_vec(:,1) = inputs.H0;
    H_o_vec(:,1) = inputs.H0*s.eta_o*rho(1)/(s.eta_o*s.old_share+s.eta_y*(1-s.old_share));
    H_y_vec(:,1) = inputs.H0-H_o_vec(:,1);
end
if isfield(inputs,'D0')
    D_vec(:,1) = inputs.D0;
    D_o_vec(:,1) = inputs.D0*varsigma(1);
    D_y_vec(:,1) = inputs.D0-D_o_vec(:,1);
end

Rt_vec = zeros(N,T); 
Rt = zeros(T,1); It = Rt; St = Rt; Iobst = Rt; Iunobst = Rt; Iobsot = Rt; Iobsyt = Rt;
Ht = Rt; Dt = Rt; Ft = Rt; %Iasympt = Rt; Isympt = Rt;
idx = ones(N,1);

% model
% S(t+1) = S(t)-R(t)*gamma*S(t)*I(t)/pop_size;
% I(t+1) = I(t)+R(t)*gamma*S(t)*I(t)/pop_size-gamma*I(t);
% z(t) = R(t)*S(t)/pop_size*<gamma,I(t)>;
for t = 1:T
    Rt_vec(:,t) = pop_size.*z(t)./S_vec(:,t).*T_si_vec./(alpha*I_unobs_vec(:,t)+I_obs_vec(:,t));
    S_vec(:,t+1) = S_vec(:,t)-z(t);
    I_unobs_vec(:,t+1) = I_unobs_vec(:,t).*(1-1./T_si_vec)+z_unobs(t);
%     I_asympt_vec(:,t+1) = I_asympt_vec(:,t).*(1-1./T_si_vec)+sigma(t).*z_obs(t); 
%     I_sympt_vec(:,t+1) = I_sympt_vec(:,t).*(1-(1-lambda)./T_si_vec-lambda./T_hosp(t))+(1-sigma(t)).*z_obs(t);    
    I_obs_y_vec(:,t+1) = I_obs_y_vec(:,t).*(1-alpha_ihy-alpha_iry)+z_obs_y(t);
    I_obs_o_vec(:,t+1) = I_obs_o_vec(:,t).*(1-alpha_iho-alpha_iro)+z_obs_o(t);
    I_obs_vec(:,t+1) = I_obs_y_vec(:,t+1)+I_obs_o_vec(:,t+1);
    H_y_vec(:,t+1) = H_y_vec(:,t).*(1-alpha_hdy-alpha_hry)+alpha_ihy.*I_obs_y_vec(:,t);
    H_o_vec(:,t+1) = H_o_vec(:,t).*(1-alpha_hdo-alpha_hro)+alpha_iho.*I_obs_o_vec(:,t);
    H_vec(:,t+1) = H_y_vec(:,t+1)+H_o_vec(:,t+1);
%     H_vec(:,t+1) = H_vec(:,t).*(1-eta_hr/T_rec-(1-eta_hr)/T_death)+lambda.*I_sympt_vec(:,t)/T_hosp(t);
%     D_vec(:,t+1) = D_vec(:,t)+(1-eta_hr)./T_death*H_vec(:,t);
    D_y_vec(:,t+1) = D_y_vec(:,t)+alpha_hdy.*H_y_vec(:,t);
    D_o_vec(:,t+1) = D_o_vec(:,t)+alpha_hdo.*H_o_vec(:,t);
    D_vec(:,t+1) = D_y_vec(:,t+1)+D_o_vec(:,t+1);
    F_vec(:,t+1) = F_vec(:,t)+alpha_hry.*H_y_vec(:,t)+alpha_hro.*H_o_vec(:,t)+...
        alpha_iry.*I_obs_y_vec(:,t)+alpha_iro.*I_obs_o_vec(:,t)+1./T_si_vec.*I_unobs_vec(:,t); 
    I_vec(:,t) = I_obs_vec(:,t)+I_unobs_vec(:,t);
    idx = idx & I_obs_vec(:,t+1)>=0 & I_unobs_vec(:,t+1)>=0 & H_vec(:,t+1)>=0;
end
idx = find(idx>0);
Rt_vec = Rt_vec(idx,:);
S_vec = S_vec(idx,:);
I_vec = I_vec(idx,:);
I_unobs_vec = I_unobs_vec(idx,:);
I_obs_vec = I_obs_vec(idx,:);
I_obs_o_vec = I_obs_o_vec(idx,:);
I_obs_y_vec = I_obs_y_vec(idx,:);
% I_asympt_vec = I_asympt_vec(idx,:);
% I_sympt_vec = I_sympt_vec(idx,:);
H_vec = H_vec(idx,:);
H_o_vec = H_o_vec(idx,:);
H_y_vec = H_y_vec(idx,:);
D_vec = D_vec(idx,:);
D_o_vec = D_o_vec(idx,:);
D_y_vec = D_y_vec(idx,:);
F_vec = F_vec(idx,:);

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

for t = 1:T
    Rt(t) = mean(Rt_vec(:,t));
    It(t) = mean(I_vec(:,t));
    Iobst(t) = mean(I_obs_vec(:,t));
    Iobsot(t) = mean(I_obs_o_vec(:,t));
    Iobsyt(t) = mean(I_obs_y_vec(:,t));
    Iunobst(t) = mean(I_unobs_vec(:,t));
%     Iasympt(t) = mean(I_asympt_vec(:,t));
%     Isympt(t) = mean(I_sympt_vec(:,t));
    St(t) = mean(S_vec(:,t));
    Ht(t) = mean(H_vec(:,t));
    Hyt(t) = mean(H_y_vec(:,t));
    Hot(t) = mean(H_o_vec(:,t));
    Dt(t) = mean(D_vec(:,t));
    Dot(t) = mean(D_o_vec(:,t));
    Dyt(t) = mean(D_y_vec(:,t));
    Ft(t) = mean(F_vec(:,t));
end

res.It = It;
res.Iot = Iobst;
res.Ioot = Iobsot;
res.Ioyt = Iobsyt;
res.Iut = Iunobst;
% res.Iat = Iasympt;
% res.Ist = Isympt;
res.Ht = Ht;
res.Hot = Hot;
res.Hyt = Hyt;
res.Dt = Dt;
res.Dot = Dot;
res.Dyt = Dyt;
res.St = St;
res.Ft = Ft;

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