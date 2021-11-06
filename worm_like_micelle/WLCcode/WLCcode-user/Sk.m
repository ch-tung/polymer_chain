function y = Sk(q,L,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
p1 = 4.12;
p2 = 4.42;
p1short = 5.36;
p2short = 5.62;
q0 = 3.1;
%qconnect = q0/b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1:length(q);

% q0short(i) = max(1.9/sqrt(Rgsquareshort(q(i),L,b)),3);
q0short(i) = max(1.9*b/sqrt(Rgsquareshort(q(i),L,b)),3);

% Modified by Guan-Rong Huang Oct. 16 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if L/b > 10
    C = 3.06/(L/b)^0.44;
else
    C = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if L/b > 10
%     epsilon = 0.176;
% else
%     epsilon = 0.170;
% end
% Modified by Guan-Rong Huang Oct. 16 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=zeros(length(q),1); % S is a 172*1 array i.e. length 172 vector
for i=1:length(q)
   if L > 4*b % Longer Chains
      if q(i)*b <= 3.1
                         
              % Modified by Yun on Oct. 15,
              Sexvmodify = Sexvnew(q, L, b);
              
              S(i) = Sexvmodify(i) + C * (4/15 + 7./(15*u(q(i),L,b)) -...
                     (11/15 + 7./(15*u(q(i),L,b))).*exp(-u(q(i),L,b)))*(b/L);
              % End of modification
                 %end
         
     else %q(i)*b > 3.1
            S(i) = a1long(q(i), L, b, p1, p2, q0)./((q(i)*b).^(p1)) + ...
                a2long(q(i), L, b, p1, p2, q0)./((q(i)*b).^(p2)) +...
                pi/(q(i)*L);
     end 
   else %L <= 4*b Shorter Chains
%      if q(i)*b <= max(1.9/sqrt(Rgsquareshort(q(i),L,b)),3)
     if q(i)*b <= max(1.9*b/sqrt(Rgsquareshort(q(i),L,b)),3) % Modified by Guan-Rong Huang Oct. 16 2019
          if q(i)*b<=0.01
              S(i) = 1 - Rgsquareshort(q(i),L,b).*(q(i).^(2))./3;
          else
           S(i) = Sdebye1(q(i),L,b);
          end
      else %q(i)*b > max(1.9/sqrt(Rgsquareshort(q(i),L,b)),3)
            S(i) = a1short(q(i),L,b,p1short,p2short,q0short(i))/((q(i)*b).^(p1short)) + ...
                a2short(q(i),L,b,p1short,p2short,q0short(i))/((q(i)*b).^(p2short)) + ...
                pi/(q(i)*L);
     end
   end
end
y = S;