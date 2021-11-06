function y = Sexv(q,L,b)
C1=1.22;
C2=0.4288;
C3=-1.651;
miu = 0.585;

y = (1 - w(q.*sqrt(Rgsquare(q,L,b)))).*Sdebye(q,L,b) + ...
    w(q.*sqrt(Rgsquare(q,L,b))).*(C1*(q.*sqrt(Rgsquare(q,L,b))).^(-1/miu) + ...
    C2*(q.*sqrt(Rgsquare(q,L,b))).^(-2/miu) +...
    C3*(q.*sqrt(Rgsquare(q,L,b))).^(-3/miu));