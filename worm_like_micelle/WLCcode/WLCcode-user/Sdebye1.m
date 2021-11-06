function y = Sdebye1(q,L,b)

y = 2*(exp(-u1(q,L,b)) + u1(q,L,b) -1)./((u1(q,L,b)).^2); 

% Checked by Guan-Rong Huang Oct. 18 2019