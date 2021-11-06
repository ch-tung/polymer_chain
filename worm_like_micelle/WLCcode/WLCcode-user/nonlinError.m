function error=nonlinError(x,jacobian,sigma)
%function error=nonlinError(x,jacobian,sigma)
%x: 
%jacobian: jacobian matrix returned by lsqnonlin or some other matlab
%fitting functions
%sigma: variance

if nargin==2
    sigma=ones(size(jacobian(:,1)));
end

N=length(x);
if N==0
    error=zeros(size(x));
    return
end

J=length(jacobian(:,1));

if length(sigma) ~= J
    printf('\sigma array is not right dimension');
    error=zeros(size(x));
    return
end


for i=1:N
    for j=1:N
        temp=sum(jacobian(:,i).*jacobian(:,j) ./ sigma);
        M(i,j)=temp(1);
    end
end

MInv = inv(M);

for i=1:N
    error(i) = sqrt(MInv(i,i));
end

return