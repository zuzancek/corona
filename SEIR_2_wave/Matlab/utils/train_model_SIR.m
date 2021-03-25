function [p]=train_model_SIR(s,data,dateFrom,dateTo)

%% scheme
% Sy' = Sy-Fy,          Fy = R/T_si*Sy*Z
% So' = So-Fo,          Fo = R/T_si*So*mu*Z
% Uy' = Uy+(1-zeta_y)*Fy-Uy/T_si;
% Uo' = Uo+(1-zeta_o)*Fo-Uo/T_si;
% Oy' = Oy+zeta_y*Fy-Oy/T_si;
% Oo' = Oo+zeta_o*Fo-Oo/T_si;
% Z = [(alpha_o*Oo+mu*Uo)/No+(alpha_y*Oy+Uy)/Ny]
%
% data
% Xobs - inflow of observed new cases
% rho - share of 65+ in new cases
% R - reproduction number (effective)

%% initialization
s.obs_ratio = 1/5;
T = dateTo-dateFrom+1;
s.sim_num = s.pop_size;
N_o = ceil(s.sim_num.*s.dep_ratio_65);
N_y = s.sim_num-N_o;
N = N_o+N_y;

method_data = s.smoothing_method_data;

% transition times
T_si = s.SI.mean+1.5; gamma = 1/T_si;
k_T = 25; sd = s.SI.std;

w_T = create_weights_time(k_T,T_si,sd); 
T_shift = max(data.from-dateFrom,k_T);

% inputs
X_obs = smooth_series(data.X_obs);
AC = method_data(data.AC);
TC = method_data(data.TC);
R_ext = resize(smooth_series(data.Rt),dateFrom:dateTo);
Rt = double(R_ext);
x_obs = double(X_obs);

% % new cases, infectious (observed/unoberved)
% X_o_obs = X_obs.*rho;                           x_o_obs = double(X_o_obs);
% X_y_obs = X_obs-X_o_obs;                        x_y_obs = double(X_y_obs);
% X_y_unobs = X_y_obs./s.obs_ratio;               x_y_unobs = double(X_y_unobs);
% X_o_unobs = s.alpha_s_o*X_o_obs./s.obs_ratio;   x_o_unobs = double(X_o_unobs);
% x_o = x_o_obs+x_o_unobs; x_o = [x_o;x_o(end);];
% x_y = x_y_obs+x_y_unobs; x_y = [x_y;x_y(end);];
% x_obs = x_o_obs+x_y_obs;
% 
% U_o = zeros(T,1); O_o = U_o; 
% U_y = zeros(T,1); O_y = U_y; 
% sigma_o = s.obs_ratio*(1+0*double(rho)*(N_y+N_o)/N_o);
% sigma_y = s.obs_ratio*(1+0*(1-double(rho))*(N_y+N_o)/N_y);
% O_o(1) = AC(dateFrom)*rho(dateFrom);            O_y(1) = AC(dateFrom)-O_o(1);
% U_o(1) = O_o(1).*s.alpha_s_o./s.obs_ratio;      U_y(1) = O_y(1)./s.obs_ratio;
% S_y = zeros(T,1);           S_y(1) = N_y-TC(dateFrom)*(1-rho(dateFrom))./s.obs_ratio;         
% S_o = zeros(T,1);           S_o(1) = N_o-TC(dateFrom)*(rho(dateFrom))./s.obs_ratio;           
% Z = S_o; Z(1) = 0;
% rt = ones(T,1); zeta = rt; 
% O = O_o; O(1) = AC(dateFrom);
% U = O; U(1) = O(1);
% alpha_o = 1;
% F = U;
% S = S_o+S_y;
% delta = .5;x_unobs = x_obs;

%% calculation
I(1) = AC(dateFrom)/s.obs_ratio;
S(1) = N-I(1);

for t=2:T_shift
    Z(t) = I(t-1)*gamma/N;   
    F(t) = Rt(t)*Z(t)*S(t-1);
    S(t) = S(t-1)-F(t);
    I(t) = I(t-1).*(1-gamma)+F(t);
end
for t=1+T_shift:T    
    Z(t) = I(t-1)*gamma/N;   
    F(t) = Rt(t)*Z(t)*S(t-1);
    S(t) = S(t-1)-F(t);
    I(t) = I(t-1)+F(t)-gamma*I(t-1);%dot(w_T(1:end-1),I(t-k_T:t-1));
end

% %% data storage
% p.S_o = tseries(dateFrom:dateTo,S_o);    p.S_y = tseries(dateFrom:dateTo,S_y);    p.S = p.S_o+p.S_y;
% p.O_o = tseries(dateFrom:dateTo,O_o);    p.O_y = tseries(dateFrom:dateTo,O_y);    p.O = p.O_o+p.O_y;
% p.U_o = tseries(dateFrom:dateTo,U_o);    p.U_y = tseries(dateFrom:dateTo,U_y);    p.U = p.U_o+p.U_y;
% p.I_o = p.O_o+p.U_o;                     p.I_y = p.O_y+p.U_y;                     p.I = p.I_o+p.I_y;
% p.sigma_o = tseries(dateFrom:dateTo,sigma_o);
% p.sigma_y = tseries(dateFrom:dateTo,sigma_y);
% p.Rt = tseries(dateFrom:dateTo,rt);
% p.zeta = tseries(dateFrom:dateTo,zeta);
% p.X_o_obs = X_o_obs;       p.X_y_obs = X_y_obs;         p.X_obs = p.X_o_obs+p.X_y_obs;
% p.X_o_unobs = X_o_unobs;   p.X_y_unobs = X_y_unobs;     p.X_unobs = p.X_o_unobs+p.X_y_unobs;
% p.X_o = p.X_o_obs+p.X_o_unobs;      p.X_y = p.X_y_obs+p.X_y_unobs;  p.X = p.X_o+p.X_y;
% 
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
plot(method_data(rt));

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

    function [pdf_x,weights] = create_weights_time(pnts_num,mean_x,std_x)
        weights = 0:pnts_num;
        pd_obj = makedist('Gamma','a',mean_x*std_x*std_x,'b',1/(std_x*std_x));
        pdf_x = pdf(pd_obj,weights);
        pdf_x = pdf_x./sum(pdf_x);
        pdf_x = pdf_x(end:-1:1);
    end

end
