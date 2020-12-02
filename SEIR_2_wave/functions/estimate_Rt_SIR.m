function [Rt,q_mat,res,Rt_last,Rt_dist,Rt_rnd] = estimate_Rt_SIR(inputs,s,do_quant,do_weight,do_dist)


% structure of inputs:
% I0: initial number of observed infectious
% z: daily data of inflow of newly observed infections

% initialization
obs_ratio = inputs.obs_ratio;
if isempty(obs_ratio)
    obs_ratio = s.obs_ratio+0*inputs.z;
end
z_obs = inputs.z;
z = z_obs./obs_ratio;
z_unobs = z-z_obs;

I0 = inputs.I0/obs_ratio(1);
N = s.sim_num;
pop_size = s.pop_size;
T = length(z);
T_si_vec = get_rv(s.SI);
alpha = ((s.T_inf_obs.mean-s.T_inf_obs0.mean)+s.T_inf_obs0.mean/s.case_isolation_effect)/s.T_inf_unobs.mean;
% set initial values
S_vec = zeros(N,T); S_vec(:,1) = pop_size-I0;
I_vec = zeros(N,T); I_vec(:,1) = I0;
I_obs_vec = zeros(N,T); I_obs_vec(:,1) = I0*obs_ratio(1);
I_unobs_vec = zeros(N,T); I_unobs_vec(:,1) = I0-I_obs_vec(:,1);
Rt_vec = zeros(N,T); 
Rt = zeros(T,1); It = Rt; St = Rt; Iobst = Rt; Iunobst = Rt;
idx = ones(N,1);

% model
% S(t+1) = S(t)-R(t)*gamma*S(t)*I(t)/pop_size;
% I(t+1) = I(t)+R(t)*gamma*S(t)*I(t)/pop_size-gamma*I(t);
% z(t) = R(t)*S(t)/pop_size*<gamma,I(t)>;
for t = 1:T
   Rt_vec(:,t) = pop_size.*z(t)./S_vec(:,t).*T_si_vec./(I_unobs_vec(:,t)+alpha*I_obs_vec(:,t));
   S_vec(:,t+1) = S_vec(:,t)-z(t);
   I_unobs_vec(:,t+1) = I_unobs_vec(:,t).*(1-1./T_si_vec)+z_unobs(t);
   I_obs_vec(:,t+1) = I_obs_vec(:,t).*(1-1./T_si_vec)+z_obs(t);
   I_vec(:,t) = I_obs_vec(:,t)+I_unobs_vec(:,t);
   idx = idx & I_unobs_vec(:,t+1)>0 & I_obs_vec(:,t+1)>0;
end
idx = find(idx>0);
Rt_vec = Rt_vec(idx,:);
S_vec = S_vec(idx,:);
I_vec = I_vec(idx,:);
I_unobs_vec = I_unobs_vec(idx,:);
I_obs_vec = I_obs_vec(idx,:);

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
    Iunobst(t) = mean(I_unobs_vec(:,t));
    St(t) = mean(S_vec(:,t));
end

res.It = It;
res.Iot = Iobst;
res.Iut = Iunobst;
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

end