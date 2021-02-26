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
date0 = min(find(resize(c,dateFrom:dateTo)>0)); %#ok<MXFND>
c(dateFrom:dateTo) = c(date0)./db_h.H(date0).*db_h.H;
db_h.C_raw = c;
db_h.C = mov_median_adj(c);
db_h.C_smooth = smooth_series(db_h.C);

% Ventilations
v = s.Ventilation;
v(dateFrom:dateTo) = v(date0)./db_h.H(date0).*db_h.H;
db_h.V_raw = v;
db_h.V = mov_median_adj(db_h.V_raw);
db_h.V_smooth = smooth_series(db_h.V);

% ECMO
e = s.OAIM;
e(dateFrom:dateTo) = e(date0)./db_h.H(date0).*db_h.H;
db_h.E_raw = e;
db_h.E = mov_median_adj(db_h.E_raw);
db_h.E_smooth = smooth_series(db_h.E);

% Serious cases (total)
db_h.S_raw = db_h.C_raw+db_h.V_raw+db_h.E_raw;
db_h.S = db_h.C+db_h.V+db_h.E;
db_h.S_smooth = smooth_series(db_h.S);

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

% ratios (severe + critical cases, mild cases vs total hospitalizations)
db_h.S_H_rate = smooth_series((db_h.C+db_h.V+db_h.E)/db_h.H);
db_h.M_H_rate = 1-db_h.S_H_rate;

% deaths
db_d.death_old_ratio_raw = a.TotalDeathRatioOld;
db_d.death_old_ratio = mov_median_adj(db_d.death_old_ratio_raw);
db_d.death_old_ratio_smooth = smooth_series(db_d.death_old_ratio);
db_d.total_raw = db_h.D_raw;        db_d.total = db_h.D;        db_d.total_smooth = db_h.D_smooth;
db_d.on_raw = db_h.D_on_raw;        db_d.on = db_h.D_on;        db_d.on_smooth = db_h.D_on_smooth;
db_d.with_raw = db_h.D_with_raw;    db_d.with = db_h.D_with;    db_d.with_smooth = db_h.D_with_smooth;
db_d.delta = mov_median_adj(db_h.D_on./db_h.D);

end
