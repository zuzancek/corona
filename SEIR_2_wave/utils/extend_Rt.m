function [Rt_new] = extend_Rt(Rt_orig,n)
[~,T] = size(Rt_orig);
Rt_ext = repmat(Rt_orig(:,T),1,n);
Rt_new = [Rt_orig,Rt_ext];
end

