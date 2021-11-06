function y = Si(q,L,b)
   
   y = zeros(size(q));
   y = sinint(L*q);

% J. S. Pedersen et al., Macromolecules 1996, 29, 7602
% Modified by Guan-Rong Huang Oct. 16 2019   
%    for i=1:length(q)
%    y(i) = quadl('sin(x)./(x)',10e-8,q(i)*L);
%    end