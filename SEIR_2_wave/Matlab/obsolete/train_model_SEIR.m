function [p]=train_model_SEIR(s,data,dateFrom,dateTo)

%% scheme
% Sy' = Sy-Fy,          Fy = R/Tinf*Sy*Z
% So' = So-Fo,          Fo = R/Tinf*So*mu*Z
% Ey' = Ey+Fy-Ey/Tlat
% Eo' = Eo+Fo-Eo/Tlat
% Uy' = Uy+(1-zeta_y)*Ey/Tlat-Uy/Tinf;
% Uo' = Uo+(1-zeta_o)*Eo/Tlat-Uo/Tinf;
% Oy' = Oy+zeta_y*E_y/Tlat-Oy/Tinf;
% Oo' = Oo+zeta_o*E_o/Tlat-Oo/Tinf;
% Z = [(alpha_o*Oo+mu*Uo)/No+(alpha_y*Oy+Uy)/Ny]
%
% data
% Xobs - inflow of observed new cases
% rho - share of 65+ in new cases
% R - reproduction number (effective)

%% initialization
s.obs_ratio = 1/3;
T = dateTo-dateFrom+1;
s.sim_num = s.pop_size;
N_o = ceil(s.sim_num.*s.dep_ratio_65);
N_y = s.sim_num-N_o;
alpha_o = s.alpha_i_o;
alpha_y = s.alpha_i_y;
alpha_s = s.alpha_s_o;

method_data = s.smoothing_method_data;

% transition times
sd = s.SI.std;
T_lat_y = s.T_lat.mean;
T_inf_y = s.T_inf.mean; gamma_y = 1/T_inf_y;
T_lat_o = s.T_lat.mean;
T_inf_o = s.T_inf.mean+2; gamma_o = 1/T_inf_o;
inf_share = s.pdf_inf_share;
inf_profile_o = s.inf_prof.fnc_s;
inf_profile_u = s.inf_prof.fnc_a;
ip_o = inf_share; %.*(1+0*inf_profile_o); 
ip_o=ip_o(end:-1:1);
ip_u = inf_share; %.*(1+0*inf_profile_u); 
ip_u=ip_u(end:-1:1);
k_inf = s.k_inf;

% weights
k_T_inf = 20; 
w_T_inf_y = create_weights_time(k_T_inf,T_inf_y,sd); 
w_T_inf_o = create_weights_time(k_T_inf,T_inf_o,sd);
k_T_lat = 25;
[~,w_T_lat_y] = create_weights_time(k_T_lat,T_lat_y,sd); 
[~,w_T_lat_o] = create_weights_time(k_T_lat,T_lat_o,sd); 

T_shift = max(k_T_inf,k_T_lat);

% inputs
rho = remove_nan(data.rho,dateFrom,dateTo);
scale = s.sim_num/s.pop_size;
X_obs = scale*smooth_series(data.X_obs);
AC = scale*method_data(data.AC);
TC = scale*method_data(data.TC);
R_ext = resize(method_data(data.Rt),dateFrom:dateTo);

% new cases, infectious (observed/unoberved)
X_o_obs = X_obs.*rho;                           x_o_obs = double(X_o_obs);
X_y_obs = X_obs-X_o_obs;                        x_y_obs = double(X_y_obs);
X_y_unobs = X_y_obs./s.obs_ratio;               x_y_unobs = double(X_y_unobs);
X_o_unobs = s.alpha_s_o*X_o_obs./s.obs_ratio;   x_o_unobs = double(X_o_unobs);
x_o = x_o_obs+x_o_unobs; x_o = [x_o;x_o(end);];
x_y = x_y_obs+x_y_unobs; x_y = [x_y;x_y(end);];

U_o = zeros(T,1); O_o = U_o; E_o = O_o;
U_y = zeros(T,1); O_y = U_y; E_y = O_y;
sigma_o = s.obs_ratio*(1+0*double(rho)*(N_y+N_o)/N_o);
sigma_y = s.obs_ratio*(1+0*(1-double(rho))*(N_y+N_o)/N_y);
O_o(1:T_shift) = AC(dateFrom:dateFrom+T_shift-1).*rho(dateFrom:dateFrom+T_shift-1).*(1-1/T_inf_o);            
O_y(1:T_shift) = AC(dateFrom:dateFrom+T_shift-1).*(1-1/T_inf_o)-O_o(1:T_shift);
U_o(1:T_shift) = O_o(1:T_shift).*(1-1/T_inf_o).*s.alpha_s_o./s.obs_ratio;      
U_y(1:T_shift) = O_y(1:T_shift).*(1-1/T_inf_o)./s.obs_ratio;
E_o(1:k_T_lat) = x_o(2:k_T_lat+1)*T_lat_o;        
E_y(1:k_T_lat) = x_y(2:k_T_lat+1)*T_lat_y;
E_o = get_wa_inv(repmat(w_T_lat_o,T+1,1),x_o,E_o(1:k_T_lat+2),1,k_T_lat+1);
E_y = get_wa_inv(repmat(w_T_lat_y,T+1,1),x_y,E_y(1:k_T_lat+2),1,k_T_lat+1);

S_y = zeros(T,1);           S_y(1) = N_y-TC(dateFrom)*(1-rho(dateFrom))./s.obs_ratio;         
S_o = zeros(T,1);           S_o(1) = N_o-TC(dateFrom)*(rho(dateFrom))./s.obs_ratio;           
Z = S_o; Z(1) = 0;
rt = ones(T,1); zeta = rt; F_o = rt; F_y = rt;

%% calculation
M = (TC(dateFrom+T_shift)-AC(dateFrom+T_shift))./s.obs_ratio;
for t=2:T_shift
    O_o(t) = O_o(t-1).*(1-gamma_o)+x_o_obs(t);
    O_y(t) = O_y(t-1).*(1-gamma_y)+x_y_obs(t);
    U_o(t) = U_o(t-1).*(1-gamma_o)+x_o_unobs(t);
    U_y(t) = U_y(t-1).*(1-gamma_y)+x_y_unobs(t);
    M=M+gamma_o*(O_o(t-1)+U_o(t-1))+gamma_y*(O_y(t-1)+U_y(t-1));
end

S_o(T_shift) = N_o-O_o(T_shift)-U_o(T_shift)-E_o(T_shift)-rho(dateFrom+T_shift)*M;
S_y(T_shift) = N_y-O_y(T_shift)-U_y(T_shift)-E_y(T_shift)-(1-rho(dateFrom+T_shift))*M;
for t=1+T_shift:T
    O_o(t) = O_o(t-1)+x_o_obs(t)-O_o(t-1)/T_inf_o;%dot(w_T_inf_o(1:end-1),O_o(t-k_T_inf:t-1));
    O_y(t) = O_y(t-1)+x_y_obs(t)-O_y(t-1)/T_inf_y;%dot(w_T_inf_y(1:end-1),O_y(t-k_T_inf:t-1));
    U_o(t) = U_o(t-1)+x_o_unobs(t)-U_o(t-1)/T_inf_o;%dot(w_T_inf_o(1:end-1),U_o(t-k_T_inf:t-1));
    U_y(t) = U_y(t-1)+x_y_unobs(t)-U_y(t-1)/T_inf_o;%dot(w_T_inf_y(1:end-1),U_y(t-k_T_inf:t-1));
    % E_o(t-1) = x_o(t)*T_lat_o;  E_o(t) = x_o(t+1)*T_lat_o;
    % E_y(t-1) = x_y(t)*T_lat_y;  E_y(t) = x_y(t+1)*T_lat_y;
    F_o(t) = E_o(t)-E_o(t-1)+x_o(t);
    F_y(t) = E_y(t)-E_y(t-1)+x_y(t);
    zeta(t) = S_y(t-1)./S_o(t-1).*F_o(t)/F_y(t);
    I_o_o = dot(ip_o(1:end-1),O_o(t-k_T_inf:t-1));
    I_y_o = dot(ip_o(1:end-1),O_y(t-k_T_inf:t-1));    
    I_o_u = dot(ip_u(1:end-1),U_o(t-k_T_inf:t-1));
    I_y_u = dot(ip_u(1:end-1),U_y(t-k_T_inf:t-1));
    Z(t) = (alpha_o*O_o(t-1)+alpha_s*U_o(t-1)).*gamma_o/N_o+...
        (alpha_y*O_y(t-1)+U_y(t-1)).*gamma_y/N_y;   
    rt(t) = F_y(t)/(S_y(t-1)*Z(t));     
    S_o(t) = S_o(t-1)-F_o(t);
    S_y(t) = S_y(t-1)-F_y(t);
end

%% data storage
p.S_o = tseries(dateFrom:dateTo,S_o);    p.S_y = tseries(dateFrom:dateTo,S_y);    p.S = p.S_o+p.S_y;
p.E_o = tseries(dateFrom:dateTo,E_o(1:end-1));    p.E_y = tseries(dateFrom:dateTo,E_y(1:end-1));    p.E = p.E_o+p.E_y;
p.O_o = tseries(dateFrom:dateTo,O_o);    p.O_y = tseries(dateFrom:dateTo,O_y);    p.O = p.O_o+p.O_y;
p.U_o = tseries(dateFrom:dateTo,U_o);    p.U_y = tseries(dateFrom:dateTo,U_y);    p.U = p.U_o+p.U_y;
p.I_o = p.O_o+p.U_o;                     p.I_y = p.O_y+p.U_y;                     p.I = p.I_o+p.I_y;
p.sigma_o = tseries(dateFrom:dateTo,sigma_o);
p.sigma_y = tseries(dateFrom:dateTo,sigma_y);
p.Rt = tseries(dateFrom:dateTo,rt);
p.zeta = tseries(dateFrom:dateTo,zeta);
p.X_o_obs = X_o_obs;       p.X_y_obs = X_y_obs;         p.X_obs = p.X_o_obs+p.X_y_obs;
p.X_o_unobs = X_o_unobs;   p.X_y_unobs = X_y_unobs;     p.X_unobs = p.X_o_unobs+p.X_y_unobs;
p.X_o = p.X_o_obs+p.X_o_unobs;      p.X_y = p.X_y_obs+p.X_y_unobs;  p.X = p.X_o+p.X_y;

%% plotting
% figure;
% subplot(2,1,1)
% plot(p.S_o,'linewidth',1); hold on;
% plot(p.S_y,'linewidth',1);
% plot(p.S,'k','linewidth',2);
% grid on;
% title('Suspectible (S)');
% legend({'Old','Young','Total'});
% 
% subplot(2,1,2)
% plot(p.E_o,'linewidth',1); hold on;
% plot(p.E_y,'linewidth',1);
% plot(p.E,'k','linewidth',2);
% grid on;
% title('Exposed (E)');
% legend({'Old','Young','Total'});
% 
% figure;
% hh1=plot(p.U_o,'linewidth',2,'linestyle',':'); hold on;
% hh2=plot(p.U_y,'linewidth',2,'linestyle',':');
% plot(p.U,'linewidth',2,'color',0.5*[1 1 1],'linestyle',':');
% plot(p.O_o,'linestyle','--','linewidth',2,'color',hh1.Color);
% plot(p.O_y,'linestyle','--','linewidth',2,'color',hh2.Color);
% plot(p.O,'linewidth',2,'color',0.5*[1 1 1],'linestyle','--');
% plot(p.I_o,'b','linewidth',3);
% plot(p.I_y,'r','linewidth',3);
% plot(p.I,'k','linewidth',3);
% grid on;
% title('Infectious (I)');
% legend({'Unobserved - Old', 'Unobserved - Young', 'Unobserved - Total',...
%     'Observed - Old', 'Observed - Young', 'Observed - Total',...
%     'Old', 'Young', 'Total'});
% 
% figure;
% subplot(2,1,1); plot(X_obs);grid on;title('X_obs');

% subplot(2,1,2); 
% plot(rt,'--');hold on;
% plot(method_data(rt));
% plot(double(R_ext));
% grid on;title('R');
hold on;
plot(smooth_series(rt));

    function [x] = remove_nan(x,t0,t1)
        if isnan(x(t0))
            i0 = find(isnan(x(t0:t1)));
            if length(i0)>1
                di0 = i0(2:end)-i0(1:end-1);
                ii0 = find(di0>1,1);
                if isempty(ii0)
                    i0 = i0(end);
                else
                    i0 = i0(ii0)-1;
                end
            else
                i0=2;
            end
            x(t0) = x(t0+i0); 
        end
        if isnan(x(t1))
            i1 = find(isnan(x(t0:t1)));
            if length(i1)>1
                di1 = i1(2:end)-i1(1:end-1);
                i1 = i1(find(di1>1,1)-1);
                if isempty(i1)
                    i1 = di1(end);
                end
            else
                i1=t1-1;
            end
            x(t1) = x(i1);
        end
        x = interp(x);   
    end

    function [pdf_x,pdf_x0,weights] = create_weights_time(pnts_num,mean_x,std_x)
        weights = 0:pnts_num;
        pd_obj = makedist('Gamma','a',mean_x*std_x*std_x,'b',1/(std_x*std_x));
        pdf_x = pdf(pd_obj,weights);
        pdf_x = pdf_x./sum(pdf_x); 
        pdf_x0 = pdf_x./weights;
        pdf_x0(1) = 0;
        pdf_x = pdf_x0(end:-1:1);
    end

    function [x] = get_wa_inv(weight,zvec,x0,alpha,idxFrom)
        sz = size(weight);
        weight = weight(idxFrom:end,:);
        k = sz(2);k0 = sz(1);
        tt = length(zvec)-idxFrom+1;        
        phi = alpha./repmat((0:k-1),tt,1); phi(:,1)=0; 
        if k0==1
            W = repmat(weight(k:-1:1),tt,1);
            A = repmat(phi(k:-1:1),tt,1);
        else
            W = weight(:,(k:-1:1));
            A = phi(:,(k:-1:1));
        end        
        J = repmat(1:k,tt,1)+repmat((0:tt-1)',1,k);
        L = repmat((1:tt)',1,k);
        Weight_mat = sparse(L(:),J(:),W(:));
        Alpha_mat = sparse(L(:),J(:),A(:));
        x = zeros(tt+k-1,1);
        x(1:k-1) = x0(1:idxFrom-1);
        function [d] = solve_lineqn(xx0)
            d=(Weight_mat.*Alpha_mat)*xx0 - zvec(end-tt+1:end);
        end
        x = max(0,fsolve(@solve_lineqn,x,optimoptions('fsolve','Display','off','Algorithm','Levenberg-Marquardt')));
    end

end
