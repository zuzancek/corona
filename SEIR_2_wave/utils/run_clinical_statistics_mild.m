function [db]=run_clinical_statistics_mild(db,s)

db_total = db{1};
db_severe = db{2};
db_mild = db_total;
ratio = s.S_H_rate;

fn = fieldnames(db_total);
n = length(fn);

for i=1:n
    pdf_total = db_total.(fn{i}).pdf;
    pdf_severe = db_severe.(fn{i}).pdf;
    k = max(length(pdf_total),length(pdf_severe));
    pdf_total(end+1:k) = 0;     db_total.(fn{i}).pdf = pdf_total;       db_total.(fn{i}).time_grid = 0:k-1;
    pdf_severe(end+1:k) = 0;    db_severe.(fn{i}).pdf = pdf_severe;     db_severe.(fn{i}).time_grid = 0:k-1;
    pdf_mild = (pdf_total-ratio.*pdf_severe)./(1-ratio);
    alpha = (db_total.(fn{i}).alpha-ratio.*db_total.(fn{i}).alpha)./(1-ratio);
    db_mild.(fn{i}).type = '';
    db_mild.(fn{i}).obj = [];
    db_mild.(fn{i}).pdf = pdf_mild;
    db_mild.(fn{i}).cdf = cumsum(pdf_mild);
    db_mild.(fn{i}).diff = 0;
    db_mild.(fn{i}).time_grid = 0:k-1;
    db_mild.(fn{i}).alpha = alpha;
end

db{1} = db_total;
db{2} = db_severe;
db{3} = db_mild;

end
