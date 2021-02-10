function [db_h,db_d] = process_clinical_statistics(s,d,a,dateFrom,dateTo)

db_h = struct; 
db_d = struct;

% hospitalizations (confirmed, suspected, total)
db_h.H_c_raw = resize(s.Hospitalizations_Confirmed,dateFrom:dateTo);
db_h.H_c = mov_median_adj(db_h.H_c_raw);
db_h.H_c_smooth = smooth_series(db_h.H_c);
db_h.H_s_raw = resize(s.Hospitalizations_Suspected,dateFrom:dateTo);
db_h.H_s = mov_median_adj(db_h.H_s_raw);
db_h.H_s_smooth = smooth_series(db_h.H_s);
db_h.H_raw = db_h.H_s_raw+db_h.H_c_raw;
db_h.H = db_h.H_s+db_h.H_c;
db_h.H_smooth = db_h.H_s_smooth+db_h.H_c_smooth;

% ICU (missing data)
c = s.ICU;
c(startdate(db_h.H_smooth)+40:startdate(db_h.H_smooth)+60) = NaN;
c = interp(c,startdate(db_h.H_smooth):enddate(db_h.H_smooth));
db_h.C_raw = c;
db_h.C = mov_median_adj(c);
db_h.C_smooth = smooth_series(db_h.C);

% Ventilations
db_h.V_raw = resize(s.Ventilation,dateFrom:dateTo);
db_h.V = mov_median_adj(db_h.V_raw);
db_h.V_smooth = smooth_series(db_h.V);

% Deaths
db_h.D_on_raw = resize(d.DeathCovid,dateFrom:dateTo);
db_h.D_on = mov_median_adj(db_h.D_on_raw);
db_h.D_on_smooth = smooth_series(db_h.D_on);
db_h.D_raw = resize(d.Total,dateFrom:dateTo);
db_h.D = mov_median_adj(db_h.D_raw);
db_h.D_smooth = smooth_series(db_h.D);
db_h.D_with_raw = db_h.D_raw-db_h.D_on_raw;
db_h.D_with = db_h.D-db_h.D_on;
db_h.D_with_smooth = db_h.D_smooth-db_h.D_on_smooth;

% Admissions
db_h.A_raw = resize(s.Admission,dateFrom:dateTo);
db_h.A = mov_median_adj(db_h.A_raw);
db_h.A_smooth = smooth_series(db_h.A);

% Discharges
db_h.R_raw = resize(s.Discharge,dateFrom:dateTo);
db_h.R = mov_median_adj(db_h.R_raw);
db_h.R_smooth = smooth_series(db_h.R);

% deaths
db_d.old_ratio_raw = a.TotalDeathRatioOld;
db_d.old_ratio = mov_median_adj(db_d.old_ratio_raw);
db_d.old_ratio_smooth = smooth_series(db_d.old_ratio);
db_d.total_raw = db_h.D_raw;        db_d.total = db_h.D;        db_d.total_smooth = db_h.D_smooth;
db_d.on_raw = db_h.D_on_raw;        db_d.on = db_h.D_on;        db_d.on_smooth = db_h.D_on_smooth;
db_d.with_raw = db_h.D_with_raw;    db_d.with = db_h.D_with;    db_d.with_smooth = db_h.D_with_smooth;

end
