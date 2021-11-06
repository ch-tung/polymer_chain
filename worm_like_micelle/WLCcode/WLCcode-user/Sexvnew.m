function y = Sexvnew(q,L,b)

C1=1.22;
C2=0.4288;
C3=-1.651;
miu = 0.585;
C_star2 = zeros(size(q));
% Modified by Guan-Rong Huang Oct. 17 2019
for i = 1:(length(q)-1)
    if Sexv(q(i+1),L,b)>=Sexv(q(i),L,b);
        C_star2(i) = 0;
    else
        C_star2(i) = 1;
    end
end
C_star2(length(q)) = C_star2(length(q)-1);

% Modified by Guan-Rong Huang Oct. 17 2019
% C_star2 = [C_star2, C_star2(length(C_star2))]';% Switch from measure q to artificial q with prime. 
y = (1 - w(q.*sqrt(Rgsquare(q,L,b)))).*Sdebye(q,L,b) + ...
    C_star2.*w(q.*sqrt(Rgsquare(q,L,b))).*(C1*(q.*sqrt(Rgsquare(q,L,b))).^(-1/miu) + ...
    C2*(q.*sqrt(Rgsquare(q,L,b))).^(-2/miu) +...
    C3*(q.*sqrt(Rgsquare(q,L,b))).^(-3/miu));