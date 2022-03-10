clear
close all

%% S_q_00
data_rod = load('scatter_rod.mat');
delta_S_q = data_rod.delta_S_q;
S_q_rod = data_rod.S_q;

filename_list = {'scatter_chain_2D_0.000100_0.000000_RB.mat',...
    'scatter_chain_2D_0.000100_0.005000_RB.mat',...
    'scatter_chain_2D_0.000100_0.010000_RB.mat',...
    'scatter_chain_2D_0.000100_0.020000_RB.mat',...
    'scatter_chain_2D_0.000100_0.040000_RB.mat'};

for i = 1:length(filename_list)
    data = load(filename_list{i});
    S_q_2D = data.S_q_2D;
    qq = data.qq;
    
    S_q_2D_plot = S_q_2D(:,:,2);
    if 1
        S_q_2D_plot = S_q_2D_plot + flipud(S_q_2D_plot) + fliplr(S_q_2D_plot) + rot90(S_q_2D_plot,2);
        S_q_2D_plot = S_q_2D_plot/4;
    end
    
    qq_2D = [-fliplr(qq),0,qq]*1000;
    
    [YY, XX] = meshgrid(qq_2D, qq_2D);
    
    % pcolor
    figure()
    pcolor(XX, YY, S_q_2D_plot)    
    colormap(parula(1000))
    caxis([0.01,1])
    
    box on
    shading interp
    set(gcf,'Position',[200,100,600,600])
    set(gca,'LineWidth',2)
    set(gca,'position',[0.24    0.22   0.72    0.72])
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    xlabel('\it{Q_x}','FontSize',24,'Interpreter','tex')
    ylabel('\it{Q_z}','FontSize',24,'Interpreter','tex')
    set(gca,'FontSize',28,'FontName','Arial')
    
    
    xlim([-2e1,2e1]*2)
    ylim([-2e1,2e1]*2)
    xticks([-10:10]*1e1)
    yticks([-10:10]*1e1)
    
    saveas(gcf,[filename_list{i},'.png'])
end

