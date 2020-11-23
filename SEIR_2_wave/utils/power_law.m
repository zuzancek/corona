function [z] = power_law(y,a0,a1,k)

    z = ((a1^(1-k)-a0^(1-k)).*y+a0^(1-k)).^(1/(1-k));

end