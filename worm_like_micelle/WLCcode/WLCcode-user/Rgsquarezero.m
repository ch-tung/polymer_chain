function y = Rgsquarezero(q,L,b)

y = ( L*b/6 ) * ( 1 - 1.5*(b/L) + 1.5*(b/L)^2 - 0.75*(b/L)^3*( 1 - exp(-2*(L/b) ) ) );

% J. S. Pedersen et al., Macromolecules 1996, 29, 7602
% Checked by Guan-Rong Huang Oct. 16 2019