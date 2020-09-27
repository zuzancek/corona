function [r] = get_R0(T,N,g1,g2,w_vec)
% g1 = reasonable/rational behaviour
% g2 = risky/irrational behaviour
r = zeros(N,T);
for t=1:T
    w = w_vec(t);
    n1 = round(w*N);
    r1 = gamrnd(g1.r0_scale*g1.mean*g1.alpha,1/g1.r0_scale,n1);
    r(1:n1,t) = r1;
    n2 = N-n1;
    r2 = power_law(rand(n2),g2.a0*g2.alpha,g2.a1*g2.alpha,g2.scale);
    r(n1+1:end,t) = r2;
end
