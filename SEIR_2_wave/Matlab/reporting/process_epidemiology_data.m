function [cases_data,test_data]=process_epidemiology_data(x,y,dateFrom,dateTo,s)

dI_inflow_ag = y.AgPosit;
dI_inflow_pcr = resize(x.NewCases,dateFrom:dateTo);
dI_inflow = dI_inflow_pcr+dI_inflow_ag;

I0 = x.TotalCases(dateFrom)/s.obs_ratio;

test_data = struct;
test_data.ptr_pcr_raw = x.NewCases./x.Tests;   test_data.ptr_pcr = mov_median_adj(test_data.ptr_pcr_raw);   test_data.ptr_pcr_smooth = smooth_series(test_data.ptr_pcr);
test_data.ptr_ag_raw = y.AgPosit./y.AgTests;   test_data.ptr_ag = mov_median_adj(test_data.ptr_ag_raw);     test_data.ptr_ag_smooth = smooth_series(test_data.ptr_ag);
test_data.tests_pcr_raw = x.Tests;             test_data.tests_pcr = mov_median_adj(x.Tests);               test_data.tests_pcr_smooth = smooth_series(test_data.tests_pcr);
d1=(startdate(test_data.ptr_ag_raw));
nanidx = find(isnan(double(test_data.ptr_ag_raw)));
test_data.tests_ag_raw = resize(y.AgTests,d1:dateTo);  
test_data.tests_ag_raw(d1+nanidx) = 0;
test_data.tests_ag = mov_median_adj(test_data.tests_ag_raw);              
test_data.tests_ag_smooth = smooth_series(test_data.tests_ag);

cases_data = struct;
cases_data.cases_ag = dI_inflow_ag;   cases_data.cases_ag_mm = mov_median_adj(dI_inflow_ag);   cases_data.cases_ag_smooth = smooth_series(cases_data.cases_ag_mm);
cases_data.cases_pcr = dI_inflow_pcr; cases_data.cases_pcr_mm = mov_median_adj(dI_inflow_pcr); cases_data.cases_pcr_smooth = smooth_series(cases_data.cases_pcr_mm);
cases_data.cases_total = dI_inflow;   cases_data.cases_total_mm = mov_median_adj(dI_inflow);   cases_data.cases_total_smooth = smooth_series(cases_data.cases_total_mm);
cases_data.I0 = I0;

end