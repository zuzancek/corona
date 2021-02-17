function [db] = process_clinical_statistics(db,s)

db_total = db{1};
db_severe = db{2};
db_s = db_total;
ratio = s.S_H_rate;

fn = fieldnames(db_total); % P_H_x, P_R_x, P_R_x ; x in {Y,O}
n = length(fn);

for i=1:n
    pdf_h = db_total.(fn{i}).pdf;       cdf_h = db_total.(fn{i}).cdf;
    pdf_s = db_severe.(fn{i}).pdf;     cdf_s = db_severe.(fn{i}).cdf;
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
    end
    
end

db{1} = db_total;
db{2} = db_m;
db{3} = db_s;

end
