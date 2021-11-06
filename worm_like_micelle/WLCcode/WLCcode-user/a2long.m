function y = a2long(q, L, b, p1, p2, q0);

if L/b > 10
    C = 3.06/(L/b)^0.44;
else
    C = 1;
end

C1 = 1.22; 
C2 = 0.4288;
C3 = -1.651;
C4 = 1.523;
C5 = 0.1477;
miu = 0.585;
m1 = -1./miu;
m2 = -2./miu;
m3 = -3./miu;
q02 = q0.^2;
q03 = q0.^3;
q04 = q0.^4;
q05 = q0.^5;
b2 = b.^2;
b3 = b.^3;
b4 = b.^4;
b5 = b.^5;
Rgs = Rgsquare(q,L,b);
z0 = (q02.*Rgs)./b2;
z1 = q02.*Rgs;
z = exp(-z0);
Rg = sqrt(Rgs);
q0Rg = Rg.*q0;
qRb = q0Rg./b;
qRb45 = (qRb - C4)./C5;

y = -q0.^(-p1).*(b2.*pi./(L.*q02) + (b.*C.*(-14.*b3./(15.*q03.*Rgs) + ...
    (14.*b3.*z)./(15.*q03.*Rgs) + (2.*z.*q0.*(11/15 + (7.*b2)./(15.*z1)).* ... 
    Rgs)./b))./L + (Rg.*(C3.*qRb.^m3 + C2.*qRb.^m2 + ...
    C1.*qRb.^m1).*sech(qRb45).^2)./(2.*C5) - ...
    (b4.*Rg.*(z + z0 -1).*sech(qRb45).^2)./(C5.*q04.*Rgs.^2) + ...
    (2.*b4.*(2*(1-z).*q0.*Rgs./b).*(1 - 0.5.*(1 + tanh(qRb45 ))))./(q04.* ... 
    Rgs.^2) - (8.*b5.*(z + z0 - 1).*(1 - 0.5.*(1 + tanh(qRb45))))./(q05.*Rgs.^2) + ...
    0.5.*(-(3.*C3.*Rg.*qRb.^(- 1 + m3))./miu - (2.*C2.*Rg.*qRb.^(- 1 + m2))./miu - ...
    (C1.*Rg.*qRb.^(-1 + m1))./miu).*(1 + tanh(qRb45))) - b.*p1.*q0.^(- 1 - p1).* ...
    (-b.*pi./(L.*q0) + (b.*C.*(4/15 - z.*(11/15 + 7.*b2./(15.*z1)) + (7.*b2)./(15.*q02.* Rgs)))./L + ...
    (2.*b4.*(z + z0 - 1).*(1 - 0.5.*(1 + tanh(qRb45))))./(q04.* Rgs.^2) + ...
    0.5.*(C3.*qRb.^m3 + C2.* qRb.^m2 + C1.*qRb.^m1).*(1 + tanh(qRb45)));

y = -1./(b.*(p1 - p2).*q0.^(-1 - p1 - p2)).*y; 

% Modified by Guan-Rong Huang Oct. 18 2019