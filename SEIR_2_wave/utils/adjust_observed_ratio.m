function [rho,rho_raw,tests,tests_raw,t0] = adjust_observed_ratio(ratio,tests,s,dateFrom)

% extend inputs
ratio = extend_series(ratio);
tests = extend_series(tests);

% initialize
rho_def = s.obs_ratio_tar;
total_tests = tests./ratio;

% calculation
d0 = min(find(resize(total_tests,s.env_from:enddate(total_tests))>s.tests_min)); %#ok<MXFND>
t0 = min(find(resize(ratio,d0:enddate(ratio))>s.ratio_threshold)); %#ok<MXFND>
ratio0 = ratio(t0);
tests0 = tests(t0);
tests_hyp = tests0.*ratio./ratio0;
scale = 1+0*tests;
scale0 = tests_hyp./tests;
scale(d0:enddate(tests)) = scale0(d0:enddate(tests));
rho_raw = rho_def./(0*scale+max(1,double(scale)));
rho = smooth_series(rho_raw,s.smooth_width,s.smooth_type,s.smooth_ends);
tests_raw = tests.*rho_def./rho_raw;
tests = smooth_series(tests_raw,s.smooth_width,s.smooth_type,s.smooth_ends);

    function [x]=extend_series(x)
        x0 = tseries(dateFrom:enddate(x),0);
        x0(startdate(x):enddate(x)) = x;
        x0(dateFrom:startdate(x)) = x0(startdate(x));
        x = x0;
    end
end