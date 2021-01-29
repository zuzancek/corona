function [db_h,db_d] = process_clinical_statistics(s,d,a,dateFrom,dateTo)

db_h = struct;
db_d = struct;

% total hospitalizations
db_h.H_raw = resize(s.Hospitalizations,dateFrom:dateTo);
db_h.H = mov_median(db_h.H_raw);
db_h.H_smooth = smooth_series(db_h.H);

% ICU (missing data)
c = s.ICU;
c(startdate(db_h.H_smooth)+40:startdate(db_h.H_smooth)+60) = NaN;
c = interp(c,startdate(hospit_smooth):enddate(hospit_smooth));
db_h.C_raw = c;
db_h.C = mov_median(c);
db_h.C_smooth = smooth_series(db_h.C);

% Ventilations
db_h.V_raw = resize(s.Ventilation,dateFrom:dateTo);
db_h.V = mov_median(db_h.V_raw);
db_h.V_smooth = smooth_series(db_h.V);

% Deaths
db_h.D_on_raw = resize(s.Deaths,dateFrom:dateTo);
db_h.D_on = mov_median(db_h.D_in_raw);
db_h.D_on_smooth = smooth_series(db_h.D_on);
db_h.D_raw = resize(d.Total,dateFrom:dateTo);
db_h.D = mov_median(db_h.D_raw);
db_h.D_smooth = smooth_series(db_h.D);
db_h.D_with_raw = db_h.D_raw-db_h.D_on_raw;
db_h.D_with = db_h.D-db_h.D_on;
db_h.D_with_smooth = db_h.D_smooth-db_h.D_on_smooth;

% Admissions
db_h.A_raw = resize(s.Admission,dateFrom:dateTo);
db_h.A = mov_median(db_h.A_raw);
db_h.A_smooth = smooth_series(db_h.A);

% Discharges
db_h.R_raw = resize(s.Discharge,dateFrom:dateTo);
db_h.R = mov_median(db_h.R_raw);
db_h.R_smooth = smooth_series(db_h.R);

% deaths
db_d.old_ratio_raw = a.TotalDeathRatioOld;
db_d.old_ratio = mov_median(db_d.old_ratio_raw);
db_d.old_ratio_smooth = smooth_series(db_d.old_ratio);
db_d.total_raw = db_d.D_raw;        db_d.total = db_h.D;        db_d.total_smooth = db_h.D_smooth;
db_d.on_raw = db_d.D_on_raw;        db_d.on = db_h.D_on;        db_d.on_smooth = db_h.D_on_smooth;
db_d.with_raw = db_d.D_with_raw;    db_d.with = db_h.D_with;    db_d.with_smooth = db_h.D_with_smooth;

end
