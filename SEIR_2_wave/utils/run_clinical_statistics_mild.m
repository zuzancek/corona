function [db]=run_clinical_statistics_mild(db,s)

db_total = db{1};
db_severe = db{2};
db_mild = db_total;
ratio = s.S_H_rate;

fn = fieldnames(db_total);
n = length(fn);

for i=1:n
    pdf_total = db_total.(fn{i}).pdf;       cdf_total = db_total.(fn{i}).cdf;
    pdf_severe = db_severe.(fn{i}).pdf;     cdf_severe = db_severe.(fn{i}).cdf;
    k = max(length(pdf_total),length(pdf_severe));
    pdf_total(end+1:k) = 0;     db_total.(fn{i}).pdf = pdf_total;       db_total.(fn{i}).time_grid = 0:k-1;
    pdf_severe(end+1:k) = 0;    db_severe.(fn{i}).pdf = pdf_severe;     db_severe.(fn{i}).time_grid = 0:k-1;
    cdf_total(end+1:k) = 1;     db_total.(fn{i}).cdf = cdf_total;       
    cdf_severe(end+1:k) = 1;    db_severe.(fn{i}).cdf = cdf_severe;  
    epdf_total(end+1:k) = 0;    db_total.(fn{i}).epdf = epdf_total;       
    epdf_severe(end+1:k) = 0;   db_severe.(fn{i}).epdf = epdf_severe;     
    ecdf_total(end+1:k) = 1;    db_total.(fn{i}).ecdf = ecdf_total;       
    ecdf_severe(end+1:k) = 1;   db_severe.(fn{i}).ecdf = ecdf_severe;  
    if strcmp(fn{i},'P_H_Y') || strcmp(fn{i},'P_H_O')
        cdf_mild = cdf_severe./cdf_total;
        pdf_mild = cdf_mild(2:end)-cdf_mild(1:end-1); pdf_mild(end+1) = 0;  pdf_mild = pdf_mild./sum(pdf_mild); %#ok<*AGROW>
        ecdf_mild = ecdf_severe./ecdf_total;
        epdf_mild = ecdf_mild(2:end)-ecdf_mild(1:end-1); epdf_mild(end+1) = 0;  epdf_mild = epdf_mild./sum(epdf_mild);
        alpha = db_severe.(fn{i}).alpha./db_total.(fn{i}).alpha;
    end
    m = (db_total.(fn{i}).mean-ratio.*db_total.(fn{i}).mean)./(1-ratio);
    db_mild.(fn{i}).type = '';
    db_mild.(fn{i}).obj = [];
    db_mild.(fn{i}).pdf = pdf_mild;
    db_mild.(fn{i}).cdf = cumsum(pdf_mild);
    db_mild.(fn{i}).epdf = epdf_mild;
    db_mild.(fn{i}).ecdf = cumsum(epdf_mild);
    db_mild.(fn{i}).diff = 0;
    db_mild.(fn{i}).time_grid = 0:k-1;
    db_mild.(fn{i}).alpha = alpha;
    db_mild.(fn{i}).mean = m;
end

db{1} = db_total;
db{2} = db_severe;
db{3} = db_mild;

end
