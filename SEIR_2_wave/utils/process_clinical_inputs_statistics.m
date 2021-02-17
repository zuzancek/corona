function [db] = process_clinical_inputs_statistics(db)

db_total = db{1};
db_severe = db{2};
db_s = db_total;
ac = {'Y','O'};

fn = {'P_H_Y','P_H_O','P_D_Y','P_D_O','P_R_Y','P_R_O'};
n = length(fn);

for i=1:n
    pdf_h = db_total.(fn{i}).pdf;       cdf_h = db_total.(fn{i}).cdf;
    pdf_s = db_severe.(fn{i}).pdf;      cdf_s = db_severe.(fn{i}).cdf;
    k = max(length(pdf_h),length(pdf_s));
    pdf_h(end+1:k) = 0;     db_total.(fn{i}).pdf = pdf_h;       db_total.(fn{i}).time_grid = 0:k-1;
    pdf_s(end+1:k) = 0;     db_severe.(fn{i}).pdf = pdf_s;      db_severe.(fn{i}).time_grid = 0:k-1;
    cdf_h(end+1:k) = 1;     db_total.(fn{i}).cdf = cdf_h;
    cdf_s(end+1:k) = 1;     db_severe.(fn{i}).cdf = cdf_s;
    epdf_h(end+1:k) = 0;    db_total.(fn{i}).epdf = epdf_h;
    epdf_s(end+1:k) = 0;    db_severe.(fn{i}).epdf = epdf_s;
    ecdf_h(end+1:k) = 1;    db_total.(fn{i}).ecdf = ecdf_h;
    ecdf_s(end+1:k) = 1;    db_severe.(fn{i}).ecdf = ecdf_s;
    if strcmp(fn{i},'P_H_Y') || strcmp(fn{i},'P_H_O')
        cdf_ms = cdf_s./cdf_h;
        pdf_ms = cdf_ms(2:end)-cdf_ms(1:end-1); pdf_ms(end+1) = 0;  pdf_ms = pdf_ms./sum(pdf_ms); %#ok<*AGROW>
        ecdf_ms = ecdf_s./ecdf_h;
        epdf_ms = ecdf_ms(2:end)-ecdf_ms(1:end-1); epdf_ms(end+1) = 0;  epdf_ms = epdf_ms./sum(epdf_ms);
        alpha = db_severe.(fn{i}).alpha./db_total.(fn{i}).alpha;
        m = dot(pdf_ms,0:k-1);
        db_s.(fn{i}).type = '';             
        db_s.(fn{i}).obj = [];
        db_s.(fn{i}).pdf = pdf_ms;
        db_s.(fn{i}).cdf = cumsum(pdf_ms);
        db_s.(fn{i}).epdf = epdf_ms;
        db_s.(fn{i}).ecdf = cumsum(epdf_ms);
        db_s.(fn{i}).diff = 0;
        db_s.(fn{i}).time_grid = 0:k-1;
        db_s.(fn{i}).alpha = alpha;
        db_s.(fn{i}).mean = m;
        db_m.(fn{i}) = db_total.(fn{i});
    elseif strcmp(fn{i},'P_D_Y') || strcmp(fn{i},'P_D_O')
        x = strcat('P_D_',ac{1+endsWith(fn{i},'O')});
        alpha = db_total.(x).alpha/db_s.(x).alpha;
        cdf_sd = db_total.(x).cdf/db_s.(x).cdf;
        pdf_sd = cdf_sd(2:end)-cdf_sd(1:end-1); pdf_sd(end+1) = 0;  pdf_sd = pdf_sd./sum(pdf_sd); %#ok<*AGROW>
        ecdf_sd = db_total.(x).ecdf/db_s.(x).ecdf;
        epdf_sd = ecdf_sd(2:end)-ecdf_sd(1:end-1); epdf_sd(end+1) = 0;  epdf_sd = epdf_sd./sum(epdf_sd);
        m = dot(pdf_sd,0:k-1);
        db_s.(fn{i}).type = '';             
        db_s.(fn{i}).obj = [];
        db_s.(fn{i}).pdf = pdf_sd;
        db_s.(fn{i}).cdf = cumsum(pdf_sd);
        db_s.(fn{i}).epdf = epdf_sd;
        db_s.(fn{i}).ecdf = cumsum(epdf_sd);
        db_s.(fn{i}).diff = 0;
        db_s.(fn{i}).time_grid = 0:k-1;
        db_s.(fn{i}).alpha = alpha;
        db_s.(fn{i}).mean = m;
        db_m.(fn{i}) = [];
    else % strcmp(fn{i},'P_R_Y') || strcmp(fn{i},'P_R_O')
        x = ac{1+endsWith(fn{i},'O')};
        alpha_m = 1-db_s.(strcat('P_H_',x)).alpha;
        alpha_s = 1-db_s.(strcat('P_D_',x)).alpha;
        cdf_sr = cdf_s; ecdf_sr = ecdf_s; pdf_sr = pdf_s; epdf_sr = epdf_s;
        cdf_mr = ((1-db_total.(strcat('P_D_',x)).alpha).*db_total.(strcat('P_D_',x)).cdf...
            -db_s.(strcat('P_H_',x)).alpha.*alpha_s.*cdf_sr.*db_s.(strcat('P_H_',x)).cdf)/alpha_m;        
        ecdf_mr = ((1-db_total.(strcat('P_D_',x)).alpha).*db_total.(strcat('P_D_',x)).ecdf...
            -db_s.(strcat('P_H_',x)).alpha.*alpha_s.*ecdf_sr.*db_s.(strcat('P_H_',x)).ecdf)/alpha_m;        
        pdf_mr = [cdf_mr(2:end)-cdf_mr(1:end-1);0]; pdf_mr = pdf_mr/sum(pdf_mr);
        epdf_mr = [ecdf_mr(2:end)-ecdf_mr(1:end-1);0]; epdf_mr = epdf_mr/sum(epdf_mr);   
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

end
