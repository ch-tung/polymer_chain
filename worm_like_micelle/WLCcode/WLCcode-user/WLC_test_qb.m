clear
close all

%%
figure()
box on
hold on

l_contour = 5000;

% load('../../scatter_chain_prstnc.mat')
load('scatter_chain_sfr_SA.mat')
plot(qq'.*a(1:4),S_q(:,1:4)./a(1:4)*l_contour,'--r','LineWidth',1)

load('scatter_chain_fvfr_SA.mat')
plot(qq'.*a(1:4),S_q(:,1:4)./a(1:4)*l_contour,'--b','LineWidth',1)

b = a;

for j = 1:4
    L = l_contour;
    bj = b(j);
    n_q = 64;
    q = (logspace(-4,0,n_q));
    
    S_m = Sk(q,L,bj);
    S_m_cs = Sk(q,L,bj).*Scs(q*bj*0.1)';
    
    plot(q*bj,S_m_cs*L/bj,':','LineWidth',1,'Color','#A0A0A0')
    plot(q*bj,S_m*L/bj,'-k','LineWidth',1)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')    
end

xlim([1e-2 1e2])
ylim([1e-1 1e3])

xlabel('$Qb$','FontSize',24,'Interpreter','latex')
ylabel('$S(Q)$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')
grid on