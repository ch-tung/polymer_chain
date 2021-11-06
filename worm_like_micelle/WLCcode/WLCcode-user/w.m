function y = w(x)

y = 0.5*( 1 + tanh( (x - 1.523)/0.1477 ) );

% C4 = 1.523;
% C5 = 0.1477;
% y = 0.5*( 1 + tanh( (x - C4)/C5 ) );
% J. S. Pedersen et al. Macromolecules 1996, 29, 7602
% Checked by Guan-Rong Huang Oct. 16 2019