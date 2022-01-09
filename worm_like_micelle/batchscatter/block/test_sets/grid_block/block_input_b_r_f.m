clear

syms x y z a b c

% p = [x*y+z - a,... 
%     -y^2*z^2+y^2*z - b,...
%     2*y*z^2-2*y*z+y - c];

q = [(z*(x+y) + (1-z)*(x)),...
    (z*(x+y-(z*(x+y) + (1-z)*(x)))^2 + (1-z)*(x-(z*(x+y) + (1-z)*(x)))^2),...
    (z*(x+y-(z*(x+y) + (1-z)*(x)))^3 + (1-z)*(x-(z*(x+y) + (1-z)*(x)))^3)/...
    ((z*(x+y-(z*(x+y) + (1-z)*(x)))^2 + (1-z)*(x-(z*(x+y) + (1-z)*(x)))^2)^(3/2))];

% p = [q(1) - a,... 
%     q(2) - b,...
%     q(3) - c];

% sol = solve(p,[x y z]);

% stats_a = [-2:0.05:-1];
% stats_b = [0:0.1:2].^2;
% stats_c = [-1:0.1:1];

params_a = [-2:0.1:-1];
params_b = [0:0.2:2];
params_c = [0.1:0.1:0.9];

[Ai,Bi,Ci] = meshgrid(params_a,params_b,params_c);
A = Ai;
B = Bi;
C = Ci;
params = {A(:),B(:),C(:)}; % [mean,var,skew]

b1 = 10.^(params{1});
rb = 10.^(params{2});
f = params{3};

parameters = [rb,b1,f];

% stats_a = eval(subs(q(1),{x y z},params));
% stats_b = eval(subs(q(2),{x y z},params));
% stats_c = eval(subs(q(3),{x y z},params));

%% plot parameter space
figure(1)
box on
scatter3(A(:),B(:),C(:),b1,'o','LineWidth',2)
xlabel('$\mu$','FontSize',24,'Interpreter','latex')
ylabel('$\sigma^2$','FontSize',24,'Interpreter','latex')
zlabel('$\gamma$','FontSize',24,'Interpreter','latex')
view([-60 18])
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])