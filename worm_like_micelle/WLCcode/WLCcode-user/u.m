function y = u(q,L,b)

y = Rgsquare(q,L,b).*(q.^2);

% J. S. Pedersen et al., Macromolecules 1996, 29, 7602
% Checked by Guan-Rong Huang