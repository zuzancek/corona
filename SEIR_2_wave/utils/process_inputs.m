function [y] = process_inputs(y,t0,t1)

fn = fieldnames(y);

if startdate(y.(fn{1}))>t0
    z = tseries(t0:t1,0);
else
    z = 0*y.(fn{1});
end

try
    for i=1:numel(fn)        
        t0y = startdate(y.(fn{i}));
        t1y = enddate(y.(fn{i}));
        x = z;
        x(t0y:t1y) = double(y.(fn{i})(t0y:t1y));
        y.(fn{i}) = x;
    end
catch er
    disp(er);
end

end