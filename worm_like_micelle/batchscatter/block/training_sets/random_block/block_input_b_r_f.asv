clear

syms x y z a b c

% p = [x*y+z - a,... 
%     -y^2*z^2+y^2*z - b,...
%     2*y*z^2-2*y*z+y - c];

q = [(z*(x+y) + (1-z)*(x)),...
    (z*(x+y-(z*(x+y) + (1-z)*(x)))^2 + (1-z)*(x-(z*(x+y) + (1-z)*(x)))^2),...
    (z*(x+y-(z*(x+y) + (1-z)*(x)))^3 + (1-z)*(x-(z*(x+y) + (1-z)*(x)))^3)/...
    ((z*(x+y-(z*(x+y) + (1-z)*(x)))^2 + (1-z)*(x-(z*(x+y) + (1-z)*(x)))^2)^(3/2))];

n_set = 10;
n_sample = 500;
n_total = n_set*n_sample;

rng(1)
params_a = rand(n_total,1)-2;
rng(2)
params_b = 2*rand(n_total,1);
rng(3)
params_c = rand(n_total,1)*0.8+0.1;

A = params_a;
B = params_b;
C = params_c;
params = {A(:),B(:),C(:)}; % [rb,b1,f]

b1 = 10.^(params{1});
rb = 10.^(params{2});
f = params{3};

parameters = [rb,b1,f];

stats_a = eval(subs(q(1),{x y z},params));
stats_b = eval(subs(q(2),{x y z},params));
stats_c = eval(subs(q(3),{x y z},params));

statistics = [stats_a,stats_b,stats_c];
%% plot parameter space
figure(1)
box on
scatter3(stats_a,stats_b,stats_c,[],b1,'o','LineWidth',2)
xlabel('$\mu$','FontSize',24,'Interpreter','latex')
ylabel('$\sigma^2$','FontSize',24,'Interpreter','latex')
zlabel('$\gamma$','FontSize',24,'Interpreter','latex')
view([-60 18])
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])