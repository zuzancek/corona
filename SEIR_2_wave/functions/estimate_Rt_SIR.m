function [Rt,q_mat,res,x_mat,Rt_last,Rt_dist,Rt_rnd] = estimate_Rt_SIR(inputs,s,do_quant,do_weight,do_dist)


% structure of inputs:
% I0: initial number of observed infectious
% z: daily data of inflow of newly observed infections

% initialization
z = inputs.z;
I0 = inputs.I0;
N = s.sim_num;
pop_size = s.pop_size;
T = length(z);
shape = s.SI.mean*(s.SI.std)^2; scale = 1/(s.SI.std)^2;
shape_vec = shape*ones(1*N,1);
scale_vec = scale*ones(1*N,1);
Trec_mat = reshape(gamrnd(shape_vec,scale_vec),N,1);

T_inf = s.T_inf_asymp;
shapeE = T_inf.mean*(T_inf.std)^2; scaleE = 1/(T_inf.std)^2;
shapeE_vec = shapeE*ones(1*N,1);
scaleE_vec = scaleE*ones(1*N,1);
Tinf_mat = reshape(gamrnd(shapeE_vec,scaleE_vec),N,1); %#ok<*NASGU>
T_inc = s.T_inc;
shapeI = T_inc.mean*(T_inc.std)^2; scaleI = 1/(T_inc.std)^2;
shapeI_vec = shapeI*ones(1*N,1);
scaleI_vec = scaleI*ones(1*N,1);
Tinc_mat = reshape(gamrnd(shapeI_vec,scaleI_vec),N,1);

% set initial values
S_vec = zeros(N,T); S_vec(:,1) = pop_size-I0;
I_vec = zeros(N,T); I_vec(:,1) = I0;
E_vec = zeros(N,T);
Rt_vec = zeros(N,T); 
X_vec = zeros(N,T);
dE_inflow = zeros(N,T);
dE_outflow = zeros(N,T);
Rt = zeros(T,1); It = Rt; Xt = Rt;Et = Rt;St = Rt;
idx = ones(N,1);
d = s.SI.mean+14;

% model
% S(t+1) = S(t)-R(t)*gamma*S(t)*I(t)/pop_size;
% I(t+1) = I(t)+R(t)*gamma*S(t)*I(t)/pop_size-gamma*I(t);
% z(t) = R(t)*gamma*S(t)*I(t)/pop_size;
for t = 1:T
   Trec_vec = Trec_mat;
   Rt_vec(:,t) = pop_size.*z(t).*Trec_vec./(S_vec(:,t).*I_vec(:,t));
   S_vec(:,t+1) = S_vec(:,t)-z(t);
   dE_inflow(:,t) = Trec_mat./T_inf.mean.*z(t); dE_outflow(:,t) = E_vec(:,t)./T_inc.mean;
   E_vec(:,t+1) = E_vec(:,t)+dE_inflow(:,t)-dE_outflow(:,t);
   I_vec(:,t+1) = I_vec(:,t).*(1-1./Trec_vec)+z(t);
   X_vec(:,t+1) = X_vec(:,t).*(1-1/d)+dE_inflow(:,t);
   idx = idx & I_vec(:,t+1)>0;
end
idx = find(idx>0);
Rt_vec = Rt_vec(idx,:);
S_vec = S_vec(idx,:);
I_vec = I_vec(idx,:);
X_vec = X_vec(idx,:);
E_vec = E_vec(idx,:);

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
    St(t) = mean(S_vec(:,t));
    Xt(t) = mean(X_vec(:,t));
    Et(t) = mean(E_vec(:,t));
end

res.It = It;
res.Iot = It*s.obs_ratio;
res.St = St;
res.Xt = Xt;
res.Et = Et;

if do_quant
    q_vec = s.quant;
    M = length(q_vec);
    q_mat = zeros(M,T);
    for j = 1:M
        q_mat(j,:) = quantile(Rt_vec,q_vec(j),1);
    end
    x_mat = zeros(M,T);
    for j = 1:M
        x_mat(j,:) = quantile(X_vec(:,2:end),q_vec(j),1);
    end
else
    q_mat = [];
    x_mat = [];
end

if do_dist
    Rt_dist = cell(T,1);
    Rt_rnd = zeros(N,T);
    % create time-varying estimator of Rt
    for t = 1:T
        r_min = min(Rt_vec(:,t));
        r_max = max(Rt_vec(:,t));
        num_pts = s.min_pts+(min(s.max_dif,max(s.min_dif,r_max-r_min))-s.min_dif)/(s.max_dif-s.min_dif)*(s.max_pts-s.min_pts);
        pts = linspace(r_min,r_max,num_pts);
        % [cdf_x,x] = ksdensity(Rt_vec(:,t),pts,'Support','positive','BoundaryCorrection','reflection','Function','cdf');
        [cdf_x,x] = ksdensity(Rt_vec(:,t),pts,'Function','cdf');
        idx = find(cdf_x(2:end)-cdf_x(1:end-1)<=0,1);
        if ~isempty(idx)
            x = x(1:idx); cdf_x = cdf_x(1:idx);
        end
        Rt_dist{t} = {x,cdf_x};
        cutoff = max([1e-5,cdf_x(1),1-cdf_x(end)]);
        u_r = random('uniform',0+cutoff,1-cutoff,N,1);
        try
            x_r = interp1(cdf_x,x,u_r,'pchip','extrap');
        catch err
            disp(t);
        end
        Rt_rnd(:,t) = reshape(x_r,[],1);
    end
else
    Rt_dist = [];
    Rt_rnd = [];
end

end