clear

syms x y z a b c

% p = [x*y+z - a,... 
%     -y^2*z^2+y^2*z - b,...
%     2*y*z^2-2*y*z+y - c];

p = [z*(x+y) + (1-z)*(x) - a,... 
    (z*(x+y-a)^2 + (1-z)*(x-a)^2) - b,...
    (z*(x+y-a)^3 + (1-z)*(x-a)^3)/(b^(3/2)) - c];

sol = solve(p,[x y z]);

n_set = 10;
n_sample = 100;
n_total = n_set*n_sample;

rng(1)
stats_a = rand(n_total,1)-2;
rng(2)
stats_b = rand(n_total,1);
rng(3)
stats_c = 2*rand(n_total,1)-1;

A = stats_a;
B = stats_b;
C = stats_c;
stats = {A(:),B(:),C(:)}; % [mean,var,skew]

lnb1_1 = eval(subs(collect(sol.x(1)),{a b c},stats));
lnrb_1 = eval(subs(collect(sol.y(1)),{a b c},stats));
f_1    = eval(subs(collect(sol.z(1)),{a b c},stats));

lnb1_2 = eval(subs(collect(sol.x(2)),{a b c},stats));
lnrb_2 = eval(subs(collect(sol.y(2)),{a b c},stats));
f_2    = eval(subs(collect(sol.z(2)),{a b c},stats));

lnb1 = lnb1_1;
lnrb = lnrb_1;
f    = f_1;

lnb1(lnrb_1<0) = lnb1_2(lnrb_1<0);
lnrb(lnrb_1<0) = lnrb_2(lnrb_1<0);
f(lnrb_1<0) = f_2(lnrb_1<0);

% b1 = exp(lnb1);
% rb = exp(lnrb);
b1 = 10.^(lnb1);
rb = 10.^(lnrb);

parameters = [rb,b1,f];
statistics = [A(:),B(:),C(:)];

%% plot parameter space
figure(1)
box on
scatter3(lnb1,lnrb,f,[],A(:),'o','LineWidth',2)
xlabel('$lnb_1$','FontSize',24,'Interpreter','latex')
ylabel('$lnr_b$','FontSize',24,'Interpreter','latex')
zlabel('$f$','FontSize',24,'Interpreter','latex')
view([-30 22.5])
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])