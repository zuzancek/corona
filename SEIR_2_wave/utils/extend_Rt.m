function [Rt_new] = extend_Rt(Rt_orig,n,varargin)
[~,T] = size(Rt_orig);
if isempty(varargin)
    Rt_ext = repmat(Rt_orig(:,T),1,n);
else
end
Rt_new = [Rt_orig,Rt_ext];
end

