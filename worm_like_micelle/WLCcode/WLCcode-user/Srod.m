function y = Srod(q,L,b)
x = q*L;
y = 2 * Si(q,L,b) ./ x - sinc(x/2/pi).^2; 

% J. S. Pedersen et al., Macromolecules 1996, 29, 7602
% Modified by Guan-Rong Huang Oct. 16 2019
% y = 2 * Si(q,L,b)./(q*L) - 4 * (sin(q*L/2).^(2))./((q*L).^(2));
