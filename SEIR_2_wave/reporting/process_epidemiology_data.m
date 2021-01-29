function [cases_data,test_data]=process_epidemiology_data(x,y,dateFrom,dateTo)

dI_inflow_ag = y.AgPosit;
dI_inflow_pcr = resize(x.NewCases,dateFrom:dateTo);
dI_inflow = dI_inflow_pcr+dI_inflow_ag;

I0 = x.TotalCases(t0)/s.obs_ratio;

test_data = struct;
test_data.ptr_pcr_raw = x.NewCases./x.Tests;   test_data.ptr_pcr = mov_median(test_data.ptr_pcr);   test_data.ptr_pcr_smooth = smooth_series(test_data.ptr_pcr);
test_data.ptr_ag_raw = y.AgPosit./y.AgTests;   test_data.ptr_ag = mov_median(test_data.ptr_ag);     test_data.ptr_ag_smooth = smooth_series(test_data.ptr_ag);
test_data.tests_pcr_raw = x.Tests;             test_data.tests_pcr = mov_median(x.Tests);           test_data.test_pcr_smooth = smooth_series(test_data.tests_pcr);
test_data.tests_ag_raw = y.AgTests;            test_data.tests_ag = mov_median(y.AgTests);          test_data.test_ag_smooth = smooth_series(test_data.tests_ag);

cases_data = struct;
cases_data.cases_ag = dI_inflow_ag;   cases_data.cases_ag_mm = mov_median(dI_inflow_ag);   cases_data.cases_ag_smooth = smooth_series(cases_data.cases_ag_mm);
cases_data.cases_pcr = dI_inflow_pcr; cases_data.cases_pcr_mm = mov_median(dI_inflow_pcr); cases_data.cases_pcr_smooth = smooth_series(cases_data.cases_pcr_mm);
cases_data.cases_total = dI_inflow;   cases_data.cases_total_mm = mov_median(dI_inflow);   cases_data.cases_total_smooth = smooth_series(cases_data.cases_total_mm);
cases_data.I0 = I0;

end