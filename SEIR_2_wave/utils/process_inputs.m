function [y] = process_inputs(y,t0,t1)

z = tseries(t0:t1,0);
fn = fieldnames(y);
try
    t0y = startdate(y.(fn{1}));
    t1y = enddate(y.(fn{1}));
    for i=1:numel(fn)
        x = z;
        x(t0y:t1y) = y.(fn{i});
        x(find(isnan(x))) = 0;
        y.(fn{i}) = x;
    end
catch er
    disp(er);
end

end