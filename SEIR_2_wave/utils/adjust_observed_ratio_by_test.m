function [ratio,ratio_adj] = adjust_observed_ratio_by_test(x,s,dateFrom,dateTo)

pos_test_ratio = x.NewCases./x.Tests;
pos_test_ratio_smooth = smooth_series(pos_test_ratio,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
tests = x.Tests;
tests_smooth = smooth_series(tests,s.smooth_width,...
    s.smooth_type,s.smooth_ends);
[tests_smooth_max,test_max_day] = max(tests_smooth);
test_ratio = tests_smooth_max./tests_smooth;
positive_ratio = pos_test_ratio_smooth./pos_test_ratio_smooth(test_max_day);
ratio = 1./(test_ratio.*positive_ratio);
ratio = resize(ratio,dateFrom:dateTo);
ratio_adj = s.obs_ratio.*ratio;

end