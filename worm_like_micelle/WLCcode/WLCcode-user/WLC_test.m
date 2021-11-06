clear
close all

%%
figure()
box on
hold on

% load('../../scatter_chain_prstnc.mat')
load('scatter_chain_sfr_woSA.mat')
plot(qq*5000,S_q(:,1:6),'--r','LineWidth',1)
qq = fliplr(qq);
S_q = fliplr(qq);

load('scatter_chain_sfr_SA.mat')
plot(qq*5000,S_q(:,1:6),'--b','LineWidth',1)
qq = fliplr(qq);
S_q = fliplr(qq);

l_contour = 5000;
b = a;

for j = 1:6
    L = l_contour;
    bj = b(j);
    n_q = 64;
    q = (logspace(-4,0,n_q));
    
    S_m = Sk(q,L,bj);
    S_m_cs = Sk(q,L,bj).*Scs(q*bj*0.1)';
    
    plot(q*L,S_m_cs,':','LineWidth',1,'Color','#A0A0A0')
    plot(q*L,S_m,'-k','LineWidth',1)
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')    
end

xlim([1e0 1e3])
ylim([2e-3 2e0])

xlabel('$Ql_{c}$','FontSize',24,'Interpreter','latex')
ylabel('$S(Q)$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Times New Roman')
grid on