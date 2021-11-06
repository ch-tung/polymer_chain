clear
close all

%% load data
% load('../../scatter_chain_prstnc.mat')
load('scatter_chain_sfr_woSA.mat')
qq = fliplr(qq);
S_q = flipud(S_q);
S_q_woSA = S_q;

load('scatter_chain_sfr_woSA.mat')
qq = fliplr(qq);
S_q = flipud(S_q);
S_q_SA = S_q;

l_contour = 5000;
b = a;

qq_max = 0.1;

figure()
box on
hold on
xi = zeros(length(b),2);
for i = 1:length(b)
    S_MC = S_q_woSA(:,i);
    S_th = @(x) Sk(qq(qq<qq_max),x(1),x(2));

    F = @(x) log(S_th(x))-log(S_MC(qq<qq_max)); % x = [L,b]
    x0 = [5000,1000];
    xlb = [0,0];
    xub = [inf,inf];

    options = optimoptions(@lsqnonlin);
    options.MaxIterations = 1000;
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(F,x0,[],[]);
    
    xi(i,:) = x;
    
    plot(qq,S_MC,'-r','LineWidth',1)
    plot(qq(qq<qq_max),S_th(x),'-k')
end

set(gca, 'XScale', 'log','LineWidth',1)
set(gca, 'YScale', 'log','LineWidth',1)

xlim([1e-4 1e-1])
ylim([2e-3 2e0])

xlabel('$Ql_{c}$','FontSize',24,'Interpreter','latex')
ylabel('$S(Q)$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')
grid on

%%
figure()
hold on
box on
axis equal
plot([1e-2 1e3],[1e-2 1e3],'--k')
plot(l_contour./b,l_contour./xi(:,2),'or')
set(gca, 'XScale', 'log','LineWidth',1)
set(gca, 'YScale', 'log','LineWidth',1)
xlim([1e-2 1e3])
ylim([1e-2 1e3])

xlabel('$l_{c}/l_{b}$','FontSize',24,'Interpreter','latex')
ylabel('$l_{c}/l_{b,fit}$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')

%%
figure()
hold on
box on
axis equal
plot([2e0 1e6],[2e0 1e6],'--k')
plot(b,xi(:,2),'or')
set(gca, 'XScale', 'log','LineWidth',1)
set(gca, 'YScale', 'log','LineWidth',1)
xlim([2e0 1e6])
ylim([2e0 1e6])

xlabel('$l_{b}$','FontSize',24,'Interpreter','latex')
ylabel('$l_{b,fit}$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')

%%
figure()
hold on
box on
% axis equal
plot([1e0 1e6],[l_contour l_contour],'--k')
plot(b,xi(:,1),'or')
set(gca, 'XScale', 'log','LineWidth',1)
% set(gca, 'YScale', 'log','LineWidth',1)
xlim([1e0 1e6])
ylim([0 1e4])

xlabel('$l_{b}$','FontSize',24,'Interpreter','latex')
ylabel('$l_{c,fit}$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')