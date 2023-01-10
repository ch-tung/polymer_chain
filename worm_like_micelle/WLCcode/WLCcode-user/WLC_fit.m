clear
close all

%% load data
% load('../../scatter_chain_prstnc.mat')
load('scatter_chain_sfr_woSA.mat')
qq = fliplr(qq);
S_q = flipud(S_q);
S_q_woSA = S_q;

load('scatter_chain_sfr_SA.mat')
qq = fliplr(qq);
S_q = flipud(S_q);
S_q_SA = S_q;

l_contour = 5000;
b = a;

qq_max = 10;

% color = parula(length(b));
color = [[0, 0.4470, 0.7410];[0.8500, 0.3250, 0.0980];[0.9290, 0.6940, 0.1250];[0.4940, 0.1840, 0.5560];[0.4660, 0.6740, 0.1880];[0.3010, 0.7450, 0.9330];	[0.6350, 0.0780, 0.1840]];

xi = zeros(length(b),2);
x_error = zeros(length(b),2);
for i = 1:length(b)
    S_MC = S_q_SA(:,i);
    S_th = @(x) Sk(qq(qq<qq_max),x(1),x(2));

    F = @(x) log(S_th(x))-log(S_MC(qq<qq_max)); % x = [L,b]
    x0 = [5000,100];
    xlb = [0,0];
    xub = [inf,inf];

    options = optimoptions(@lsqnonlin);
    options.MaxIterations = 1000;
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(F,x0,[],[]);
    
    xi(i,:) = x;
    
    figure(1)
    box on
    hold on
    plot(qq(1:2:end),S_MC(1:2:end),'o','LineWidth',1,'Color',color(i,:))
    plot(qq(qq<qq_max),S_th(x),'-','Color',color(i,:))
    
    figure(2)
    box on
    hold on
    plot(qq(1:2:end)*b(i),S_MC(1:2:end)*l_contour/b(i),'o','LineWidth',1,'Color',color(i,:))
    plot(qq(qq<qq_max)*b(i),S_th(x)*l_contour/b(i),'-','Color',color(i,:))
    
    x_error(i,:) = nonlinError(x,jacobian);
end

figure(1)
set(gca, 'XScale', 'log','LineWidth',1)
set(gca, 'YScale', 'log','LineWidth',1)

xlim([1e-4 1e-1])
ylim([2e-3 2e0])

xlabel('$Ql_{c}$','FontSize',24,'Interpreter','latex')
ylabel('$S(Q)$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')
% grid on

figure(2)
set(gca, 'XScale', 'log','LineWidth',1)
set(gca, 'YScale', 'log','LineWidth',1)

xlim([1e-2 1e2])
ylim([1e-1 1e3])

xlabel('$Ql_{b}$','FontSize',24,'Interpreter','latex')
ylabel('$(l_{c}/l_{b}) S(Q)$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')
% grid on

%%
figure(3)
hold on
box on
axis equal
plot([5e-1 2e3],[5e-1 2e3],'--k')
plot(l_contour./b,l_contour./xi(:,2),'or','MarkerSize',12,'LineWidth',1)
set(gca, 'XScale', 'log','LineWidth',1)
set(gca, 'YScale', 'log','LineWidth',1)
xlim([5e-1 2e3])
ylim([5e-1 2e3])

xlabel('$l_{c}/l_{b}$','FontSize',24,'Interpreter','latex')
ylabel('$l_{c}/l_{b,fit}$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')

%%
figure(4)
hold on
box on
axis equal
plot([1e0 1e4],[1e0 1e4],'--k')
plot(b,xi(:,2),'or','MarkerSize',12,'LineWidth',1)
set(gca, 'XScale', 'log','LineWidth',1)
set(gca, 'YScale', 'log','LineWidth',1)
xlim([1e0 1e4])
ylim([1e0 1e4])

xlabel('$l_{b}$','FontSize',24,'Interpreter','latex')
ylabel('$l_{b,fit}$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')

%%
figure(5)
hold on
box on
% axis equal
plot([1e-1 1e4],[l_contour l_contour],'--k')
plot(l_contour./b,xi(:,1),'or','MarkerSize',12,'LineWidth',1)
set(gca, 'XScale', 'log','LineWidth',1)
% set(gca, 'YScale', 'log','LineWidth',1)
xlim([5e-1 2e3])
ylim([0 1e4])

xticks(10.^[0:3])

xlabel('$l_{c}/l_{b}$','FontSize',24,'Interpreter','latex')
ylabel('$l_{c,fit}$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')