function [ratio_smooth,ratio_adj,new_cases] = adjust_observed_ratio_by_test(x,s,dateFrom,dateTo)

d0 = min(find(resize(x.Tests,s.env_from:enddate(x.Tests))>s.tests_min))-30; %#ok<MXFND>

pos_test_ratio = x.NewCases./x.Tests;
pos_test_ratio_smooth = smooth_series(pos_test_ratio,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
tests = x.Tests;
tests_smooth = smooth_series(tests,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
[tests_smooth_max,test_max_day] = max(tests_smooth);
test_ratio = tests_smooth_max./tests_smooth;
test_ratio(dateFrom:d0-1) = 1;
%positive_ratio = pos_test_ratio_smooth./pos_test_ratio_smooth(test_max_day);
ratio = 1./(test_ratio);

d0 = min(find(resize(x.Tests,s.env_from:enddate(x.Tests))>s.tests_min))-30; %#ok<MXFND>
ratio(dateFrom:d0-1) = 1;
%ratio = ratio.*(x.Tests>=2000 & pos_test_ratio>=0.05)+1;
ratio = resize(ratio,dateFrom:dateTo);
ratio_smooth = ratio;
%smooth_series(ratio,s.smooth_width,...
%    s.smooth_type,s.smooth_ends);
ratio_adj = s.obs_ratio.*ratio_smooth;

new_cases = resize(x.NewCases,dateFrom:dateTo);
new_cases = smooth_series(new_cases,s.smooth_width,s.smooth_type,s.smooth_ends);
new_cases = new_cases./ratio_smooth;
new_cases = smooth_series(new_cases,s.smooth_width,s.smooth_type,s.smooth_ends);

end