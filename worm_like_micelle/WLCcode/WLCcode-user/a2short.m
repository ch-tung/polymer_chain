function y = a2short(q, L, b, p1, p2, q0);

Rgs = Rgsquareshort(q,L,b);
q02 = q0.^2;
q03 = q0.^3;
b2 = b.^2;
b3 = b.^3;
z = exp(q02.*Rgs./b2);

y = 8.*b3.*L - 8.*b3.*z.*L - 2.*b3.*L.*p1 + 2.*b3.*z.*L.*p1 + ...
    4.*b.*L.*q02.*Rgs + 4.*b.*z.*L.*q02.*Rgs - 2.*b.*z.*L.*p1.*q02.*Rgs - ... 
    z.*pi.*q03.*Rgs.^2 + z.*p1.*pi.*q03.*Rgs.^2; 

y = -b./z.*q0.^(p2 - 4).*y./(L.*(p1 - p2).*Rgs.^2);

% Modified by Guan-Rong Huang Oct. 18 2019;