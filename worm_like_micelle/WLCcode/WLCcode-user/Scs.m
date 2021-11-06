function y = Scs(x)

k = find(x==0);
y = (2.*besselj(1,x)./x).^(2);
y(k) = 1;

% Modified by Guan-Rong Huang Oct. 16 2019