clear
close all

%%
figure()
box on
hold on

l_contour = 5000;


L = l_contour;
bj = 100;
n_q = 64;
q = (logspace(-4,0,n_q));

S_m = Sk(q,L,bj);
S_m_cs = Sk(q,L,bj).*Scs(q*bj*0.1)';

plot(q*bj,S_m,'-k','LineWidth',2)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')

xlim([1e-2 1e2])
ylim([2e-4 2e0])

xlabel('$kl_p$','FontSize',24,'Interpreter','latex')
ylabel('$S(k)$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,800,600])
set(gca,'FontSize',28,'FontName','Arial')
% grid on