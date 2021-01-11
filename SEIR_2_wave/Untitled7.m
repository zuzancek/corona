close all;
k = 15;
x = 0:k;
xx = 0:0.01:k;
d_std = 0.62; d_mean = 6.5;
d_shape = d_mean*d_std^2; d_scale = 1/d_std^2;
y = pdf('Gamma',x,d_shape,d_scale);
yy = pdf('Gamma',xx,d_shape,d_scale);
figure;
plot(xx,yy);hold on;
plot(x,y,'o')
s = sum(y);
y = y/s;
s = '';
for i=1:length(x)
    s = strcat(s,'%2.2f\t');
end
s = strcat(s,'\n');
fprintf(s,x);
fprintf(s,100*y);
%
T = 100;
w = y(2:end);
V = repmat(w,T,1);
J = repmat([1:k],T,1)+repmat([0:T-1]',1,k);
I = k-1+repmat([1:T]',1,k);
A0 = tril(repmat([1:k-1],k-1,1)); 
A0(A0==0) = k+1;
w(end+1) = 0;
cw = cumsum(w);
J0 = repmat(1:k-1,k-1,1);
I0 = repmat([1:k-1]',1,k-1);
V0 = w(A0(A0~=0));
A = sparse([I0(:);I(:)],[J0(:);J(:)],[V0(:);V(:)]);
A = A./sum(A,2);
B = full(A);