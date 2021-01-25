function [yvec]=eval_at(tvec,xvec,tmax)

tidx = min(floor(tvec),tmax-1);
yvec = (tvec-tidx).*xvec(tidx+1)+(tidx+1-tvec).*xvec(tidx);

end