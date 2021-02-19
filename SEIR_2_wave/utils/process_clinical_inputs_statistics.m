function [db] = process_clinical_inputs_statistics(db)

db_total = db{1};
db_severe = db{2};
db_s = db_total;
ac = {'y','o'};

fn = strcat('opt_fit_',{'h_y','h_o','d_y','d_o','r_y','r_o'});
n = length(fn);

for i=1:n
    pdf_h = db_total.(fn{i}).pdf;       cdf_h = db_total.(fn{i}).cdf;
    pdf_s = db_severe.(fn{i}).pdf;      cdf_s = db_severe.(fn{i}).cdf;
    epdf_h = db_total.(fn{i}).epdf;     ecdf_h = db_total.(fn{i}).ecdf;
    epdf_s = db_severe.(fn{i}).epdf;    ecdf_s = db_severe.(fn{i}).ecdf;
    k = max(length(pdf_h),length(pdf_s));
    pdf_h(end+1:k) = 0;     db_total.(fn{i}).pdf = pdf_h;       db_total.(fn{i}).time_grid = 0:k-1;
    pdf_s(end+1:k) = 0;     db_severe.(fn{i}).pdf = pdf_s;      db_severe.(fn{i}).time_grid = 0:k-1;
    cdf_h(end+1:k) = 1;     db_total.(fn{i}).cdf = cdf_h;
    cdf_s(end+1:k) = 1;     db_severe.(fn{i}).cdf = cdf_s;
    epdf_h(end+1:k) = 0;    db_total.(fn{i}).epdf = epdf_h;
    epdf_s(end+1:k) = 0;    db_severe.(fn{i}).epdf = epdf_s;
    ecdf_h(end+1:k) = 1;    db_total.(fn{i}).ecdf = ecdf_h;
    ecdf_s(end+1:k) = 1;    db_severe.(fn{i}).ecdf = ecdf_s;
    if endsWith(fn{i},'h_y') || endsWith(fn{i},'h_o')
        cdf_ms = adjust_cdf(min(1,max(0,cdf_s./cdf_h))); cdf_ms(isnan(cdf_ms)) = 0;
        pdf_ms = [0;cdf_ms(2:end)-cdf_ms(1:end-1)];  pdf_ms = pdf_ms./sum(pdf_ms); %#ok<*AGROW>
        ecdf_ms = ecdf_s./ecdf_h; 
        ecdf_ms = adjust_cdf(ecdf_ms);
        epdf_ms = [0;ecdf_ms(2:end)-ecdf_ms(1:end-1)];  epdf_ms = epdf_ms./sum(epdf_ms);
        alpha = db_severe.(fn{i}).alpha./db_total.(fn{i}).alpha;
        m = dot(pdf_ms,0:k-1);
        db_s.(fn{i}).type = '';             
        db_s.(fn{i}).obj = [];
        db_s.(fn{i}).pdf = pdf_ms;
        db_s.(fn{i}).cdf = cdf_ms;
        db_s.(fn{i}).epdf = epdf_ms;
        db_s.(fn{i}).ecdf = ecdf_ms;
        db_s.(fn{i}).diff = 0;
        db_s.(fn{i}).time_grid = 0:k-1;
        db_s.(fn{i}).alpha = alpha;
        db_s.(fn{i}).mean = m;
        db_m.(fn{i}) = db_total.(fn{i});
    elseif endsWith(fn{i},'d_y') || endsWith(fn{i},'d_o')
        x0 = ac{1+endsWith(fn{i},'o')};        
        xd = strcat('opt_fit_d_',x0);       
        xh = strcat('opt_fit_h_',x0);
        db_s.(xd) = db_severe.(xd);
        alpha = db_total.(xd).alpha/db_s.(xh).alpha;
        cdf_sd = min(1,max(0,db_total.(xd).cdf./db_s.(xh).cdf));  cdf_sd(isnan(cdf_sd)) = 0;
        pdf_sd = max(0,[0;cdf_sd(2:end)-cdf_sd(1:end-1)]);  pdf_sd = pdf_sd./sum(pdf_sd); %#ok<*AGROW>
        ecdf_sd = min(1,max(0,db_total.(xd).ecdf./db_s.(xh).ecdf)); ecdf_sd(isnan(ecdf_sd)) = 0;
        epdf_sd = max(0,[0;ecdf_sd(2:end)-ecdf_sd(1:end-1)]);  epdf_sd = epdf_sd./sum(epdf_sd);
        m = dot(pdf_sd,0:k-1);
        db_s.(fn{i}).type = '';             
        db_s.(fn{i}).obj = [];
        db_s.(fn{i}).pdf = pdf_sd;
        db_s.(fn{i}).cdf = cdf_sd;
        db_s.(fn{i}).epdf = epdf_sd;
        db_s.(fn{i}).ecdf = ecdf_sd;
        db_s.(fn{i}).diff = 0;
        db_s.(fn{i}).time_grid = 0:k-1;
        db_s.(fn{i}).alpha = alpha;
        db_s.(fn{i}).mean = m;
        db_m.(fn{i}) = [];
    else % strcmp(fn{i},'P_R_Y') || strcmp(fn{i},'P_R_O')
        x0 = ac{1+endsWith(fn{i},'o')};        
        xd = strcat('opt_fit_d_',x0);       
        xh = strcat('opt_fit_h_',x0);      
        xr = strcat('opt_fit_r_',x0);
        db_s.(xr) = db_severe.(xr);
        alpha_m = 1-db_s.(xh).alpha;
        alpha_s = 1-db_s.(xd).alpha;
        alpha_r = db_total.(xr).alpha;
        cdf_sr = cdf_s; ecdf_sr = ecdf_s; pdf_sr = pdf_s; epdf_sr = epdf_s;
        cdf_mr = min(1,max(0,(alpha_r.*db_total.(xr).cdf...
            -db_s.(xh).alpha.*alpha_s.*cdf_sr.*db_s.(xh).cdf)/alpha_m));        
        ecdf_mr = min(1,max(0,(alpha_r.*db_total.(xr).ecdf...
            -db_s.(xh).alpha.*alpha_s.*ecdf_sr.*db_s.(xh).ecdf)/alpha_m));        
        pdf_mr = max(0,[0;cdf_mr(2:end)-cdf_mr(1:end-1)]); pdf_mr = pdf_mr/sum(pdf_mr);
        epdf_mr = max(0,[0;ecdf_mr(2:end)-ecdf_mr(1:end-1)]); epdf_mr = epdf_mr/sum(epdf_mr);   
        m_m = dot(epdf_mr,0:k-1);
        m_s = dot(epdf_sr,0:k-1);
        db_s.(fn{i}) = db_severe.(fn{i});
        db_s.(fn{i}).time_grid = 0:k-1;
        db_m.(fn{i}).time_grid = 0:k-1;
        db_s.(fn{i}).alpha = alpha_s;
        db_m.(fn{i}).alpha = alpha_m;
        db_s.(fn{i}).pdf = pdf_sr;
        db_s.(fn{i}).epdf = epdf_sr;
        db_s.(fn{i}).cdf = cdf_sr;
        db_s.(fn{i}).ecdf = ecdf_sr;
        db_m.(fn{i}).pdf = pdf_mr;
        db_m.(fn{i}).epdf = epdf_mr;
        db_m.(fn{i}).cdf = cdf_mr;
        db_m.(fn{i}).ecdf = ecdf_mr; 
        db_m.(fn{i}).type = '';             
        db_m.(fn{i}).obj = [];
        db_m.(fn{i}).diff = 0;
        db_m.(fn{i}).mean = m_m;
        db_s.(fn{i}).mean = m_s;        
    end    
end

db{1} = db_total;
db{2} = db_m;
db{3} = db_s;

    function [z]=adjust_cdf(y)
        dy = y(2:end)-y(1:end-1);
        idx = find(dy<0);
        if ~isempty(idx)
            y(2:idx(end)+1) = NaN; y(1) = 0;
            z = (reshape(interp1(find(~isnan(y)),y(find(~isnan(y))),1:length(y),'pchip'),[],1)); %#ok<FNDSB>
        else
            z = y;
        end
    end

end

