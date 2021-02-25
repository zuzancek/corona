function [Rt,q_mat,res,Rt_last,Rt_dist,Rt_rnd] = estimate_Rt_SIR(inputs,s,do_quant,do_weight,do_dist)

%% Model
% S(t+1) = S(t)-Z(t)
% Iu(t+1) = Iu(t)+Z(t)-rho(t)*Iu(t)/T_test-(1-rho(t))*Iu(t)/T_inf
% Io(t+1) = Io(t)+rho(t)*Iu(t)/T_test-lambda(t)*Io(t)/T_hosp-(1-lambda(t))*Io(t)/T_sick
%%

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
    eta_y = inputs.eta_y;
    eta_y = extend(eta_y,T-length(eta_y));
catch err %#ok<NASGU>
    eta_y = s.eta_y+0*inputs.z;
end
try
    eta_o = inputs.eta_o;
    eta_o = extend(eta_o,T-length(eta_o));
catch err %#ok<NASGU>
    eta_o = s.eta_o+0*inputs.z;
end
try
    rho = double(inputs.old_ratio);
    if length(rho)>=T
        rho = rho(end-T+1:end);
    else
        rho = extend(rho,T-length(rho));
    end
catch err %#ok<NASGU>
    rho = s.old_share+0*inputs.z;
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
alpha_ihy = s.eta_y/s.T_hosp_y_mean;        
alpha_iho = s.eta_o/s.T_hosp_o_mean;  
x.T_sick_y.mean = s.T_sick_y-s.T_test.mean; x.T_sick_y.std = s.T_sick_std;
T_sick_y_vec = get_rv(x.T_sick_y);
alpha_iry_vec = (1-eta_y)'./T_sick_y_vec;
x.T_sick_o.mean = s.T_sick_o-s.T_test.mean; x.T_sick_o.std = s.T_sick_std;
T_sick_o_vec = get_rv(x.T_sick_o);
alpha_iro_vec = (1-eta_o)'./T_sick_o_vec;
% T_hosp = s.T_hosp.mean;
alpha = 1.15; %((s.T_inf_obs.mean-s.T_inf_obs0.mean)+s.T_inf_obs0.mean/s.case_isolation_effect)/s.T_inf_unobs.mean;
% set initial values
S_vec = zeros(N,T); S_vec(:,1) = pop_size-I0;
I_vec = zeros(N,T); I_vec(:,1) = I0;
I_obs_vec = zeros(N,T); I_obs_vec(:,1) = I0*obs_ratio(1);
I_unobs_vec = zeros(N,T); I_unobs_vec(:,1) = I0-I_obs_vec(:,1);
I_obs_o_vec = zeros(N,T); I_obs_o_vec(:,1) = I0*s.old_share;
I_obs_y_vec = zeros(N,T); I_obs_y_vec(:,1) = I0-I_obs_o_vec(:,1);
% I_asympt_vec = zeros(N,T); I_asympt_vec(:,1) = I_obs_vec(:,1).*sigma(1);
% I_sympt_vec = zeros(N,T); I_sympt_vec(:,1) = I_obs_vec(:,1).*(1-sigma(1));

Rt_vec = zeros(N,T); 
Rt = zeros(T,1); It = Rt; St = Rt; Iobst = Rt; Iunobst = Rt; Iobsot = Rt; Iobsyt = Rt;
%Iasympt = Rt; Isympt = Rt;
idx = ones(N,1);

% model
% S(t+1) = S(t)-R(t)*gamma*S(t)*I(t)/pop_size;
% I(t+1) = I(t)+R(t)*gamma*S(t)*I(t)/pop_size-gamma*I(t);
% z(t) = R(t)*S(t)/pop_size*<gamma,I(t)>;
for t = 1:T
    Rt_vec(:,t) = pop_size.*z(t)./S_vec(:,t).*T_si_vec./(alpha*I_unobs_vec(:,t)+I_obs_vec(:,t));
    S_vec(:,t+1) = S_vec(:,t)-z(t);
    I_unobs_vec(:,t+1) = I_unobs_vec(:,t).*(1-1./T_si_vec)+z_unobs(t);  
    I_obs_y_vec(:,t+1) = I_obs_y_vec(:,t).*(1-alpha_ihy-alpha_iry_vec(:,t))+z_obs_y(t);
    I_obs_o_vec(:,t+1) = I_obs_o_vec(:,t).*(1-alpha_iho-alpha_iro_vec(:,t))+z_obs_o(t);
    I_obs_vec(:,t+1) = I_obs_y_vec(:,t+1)+I_obs_o_vec(:,t+1);
     
    I_vec(:,t) = I_obs_vec(:,t)+I_unobs_vec(:,t);
    idx = idx & I_obs_vec(:,t+1)>=0 & I_unobs_vec(:,t+1)>=0;
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
end

res.It = It;
res.Iot = Iobst;
res.Ioot = Iobsot;
res.Ioyt = Iobsyt;
res.Iut = Iunobst;
% res.Iat = Iasympt;
% res.Ist = Isympt;
res.St = St;

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

    function [y] = extend(x,t0)
        [xlen,xwid] = size(x);
        z = x(1,:)+zeros(xlen+t0,xwid);
        z(t0+1:end,:) = x;
        y = (z);
    end

end