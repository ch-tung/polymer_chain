function y = Sdebye(q,L,b)

y = 2*(exp(-u(q,L,b)) + u(q,L,b) -1) ./ ( (u(q,L,b)).^2 ); 

% J. S. Pedersen et al., Macromolecules 1996, 29, 7602
% Checked by Guan-Rong Huang Oct. 16 2019