clear
close all

%% S_q_00
data_rod = load('scatter_rod.mat');
delta_S_q = data_rod.delta_S_q;
S_q_rod = data_rod.S_q;

filename_list = {'scatter_chain_aniso_0.000100_0.000000_RB.mat',...
    'scatter_chain_aniso_0.000100_0.005000_RB.mat',...
    'scatter_chain_aniso_0.000100_0.010000_RB.mat',...
    'scatter_chain_aniso_0.000100_0.020000_RB.mat',...
    'scatter_chain_aniso_0.000100_0.040000_RB.mat'};

for i = 1:length(filename_list)
    data = load(filename_list{i});
    S_q_00 = data.S_q_lm(:,1);
    qq = data.qq;
    
    figure(1)
    hold on
    
    plot(qq*1000,S_q_00'.*delta_S_q)
    
    box on
    set(gcf,'Position',[200,100,600,600])
    set(gca,'LineWidth',2)
    set(gca,'position',[0.22    0.22   0.72    0.72])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('\it{QL}','FontSize',24,'Interpreter','tex')
    ylabel('{\it S}({\itQL})','FontSize',24,'Interpreter','tex')
    set(gca,'FontSize',28,'FontName','Arial')
    
    
    xlim([1e-1,1e3])
    ylim([1e-3,2e0])
    xticks([1e-1,1e0,1e1,1e2,1e3])
    yticks([1e-3,1e-2,1e-1,1e0])
end

plot(qq*1000,pi./qq/1000,'--k')

for i = 1:length(filename_list)
    data = load(filename_list{i});
    S_q_00 = data.S_q_lm(:,1);
    qq = data.qq;
    
    figure(2)
    hold on
    
    plot(qq*1000,S_q_00'.*delta_S_q.*qq*1000)
    
    box on
    set(gcf,'Position',[200,100,600,600])
    set(gca,'LineWidth',2)
    set(gca,'position',[0.22    0.22   0.72    0.72])
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('\it{QL}','FontSize',24,'Interpreter','tex')
    ylabel('{\it QLS}({\itQL})','FontSize',24,'Interpreter','tex')
    set(gca,'FontSize',28,'FontName','Arial')
    
    
    xlim([1e-1,1e3])
    ylim([1e-1,2e1])
    xticks([1e-1,1e0,1e1,1e2,1e3])
    yticks([1e-3,1e-2,1e-1,1e0,1e1])
end

plot(qq*1000,pi./qq.*qq,'--k')
